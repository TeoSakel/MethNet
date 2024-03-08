#!/usr/bin/env Rscript
# Parse inputs -----------------------------------------------------------------
doc <- "
Runs elastic net model for every gene independently.

Usage: fit_enet_models.R [options] <net> <meth> <gex> <prefix>

Options:
  --clusters <PATH>                file mapping probe_ids to cluster_ids
  --cluster <COLNAME>              which cluster scheme to use
  -l <LAMBDA>, --lambda <LAMBDA>   which lambda to use [default: lambda.1se]
  --family <FAMILY>                GLM family to use [default: gaussian]
  --norm                           Normalize gene expression before fitting.
  --gex_var=<VALUE>                filter for gene expression variance [default: 0]
  --gex_uq=<VALUE>                 filter for gene expression variance [default: 0]
  --samples <PATH>                 file of samples to use
  --genes <PATH>                   file of genes to fit
  --clinical                       use clinical information in the model
"
library(docopt)
argv <- docopt(doc)

# crash early
for (arg in c("net", "meth", "gex", "clusters", "samples", "genes")) {
    if (!is.null(argv[[arg]]) && !file.exists(argv[[arg]]))
        stop(sprintf("File %s does not exist", argv[[arg]]))
}

# required
net_file  <- argv[["net"]]
meth_file <- argv[["meth"]]
gex_file  <- argv[["gex"]]
prefix    <- argv[["prefix"]]
# optional
lambda         <- argv[["lambda"]]
family         <- argv[["family"]]
scale_gex      <- argv[["norm"]]
gex_min_mad    <- as.numeric(argv[["gex_var"]])
gex_min_uq     <- as.numeric(argv[["gex_uq"]])
cluster_file   <- argv[["clusters"]]
cluster_scheme <- argv[["cluster"]]
samples_file   <- argv[["samples"]]
genes_file     <- argv[["genes"]]
use_clinical   <- argv[["clinical"]]

# sanity checks
stopifnot(lambda %in% c("lambda.min", "lambda.1se"))
stopifnot(all(!is.na(c(gex_min_mad, gex_min_uq))))

prefix_dir <- dirname(prefix)
if (!dir.exists(prefix_dir)) {
    warning(sprintf("Creating directory %s to write results", prefix_dir))
    dir.create(prefix_dir, recursive = TRUE)
}
# fix prefix (to work with snakemake)
cancer <- sub("_.*", "", basename(prefix))
prefix <- file.path(prefix_dir, cancer)

# setup parallel (through environment inputs)
suppressPackageStartupMessages(library(parallel))
ncores <- future::availableCores()
options("mc.cores" = ncores)


# Setup environment ------------------------------------------------------------

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(glmnet)
    library(glmnetUtils)
    library(plotmo)
})

fread_wrapper <- function(x) {
    if (tools::file_ext(x) == "rds") {
        readRDS(x)
    } else if (grepl(".gz$", x)) {
        fread(cmd = paste("zcat", x))
    } else {
        fread(x)
    }
}

get_cva_attr <- function(mod, s = "lambda.1se") {
    attrs <- c("cvm", "cvsd", "cvup", "cvlo", "nzero")
    if (s %in% c("lambda.min", "lambda.1se")) {
        lam <- mod[[s]]
    } else if (is.numeric(s)) {
        lam <- mod[["lambda"]]
        dlam <- abs(lam - s)
        ix <- which.min(dlam)
        lam <- lam[ix]
        if (dlam[ix] > 1) {
            warning(paste("s is not close to any lambda value.",
                          "Better specify fixed lambda during training"))
            return(NA_real_)
        } else if (dlam[ix] > 0.001) {
            warning(sprintf("s = %.3f is not availabe as an option. Using s = %.3f instead.", s, lam))
        }
    } else {
        stop("Argument `s` must be either numeric or one of 'lambda.min', 'lambda.1se'")
    }
    # top-level attributes
    ix <- which(mod[["lambda"]] == lam)
    res <- vapply(attrs, function(x) as.numeric(mod[[x]][ix]), 0)
    # number of non-zero coefs exluding intercept:
    res["nnzero"] <- mod$glmnet.fit$dim[1L] - res["nzero"]
    # attributes inside glmnet.fit
    ix <- which(mod$glmnet.fit$lambda == lam)  # sometimes it's not the same
    res["R2"] <- mod$glmnet.fit$dev.ratio[ix]
    res[c("R2", attrs, "nnzero")]
}

get_cva_score <- function(mod, s = "lambda.1se") {
    get_cva_attr(mod, s)["cvm"]
}

# TODO: write to disk to avoid memory drain?
if (use_clinical) {
    prepare_matrix <- function(gene) {
        frml <- as.formula(sprintf(" ~ (1 + %s) * clinical", paste(net[[gene]], collapse = " + ")))
        X <- meth[, net[[gene]], drop = FALSE]
        X <- model.matrix(frml, data = cbind(data.frame("clinical" = clinical), as.data.frame(X)))
        X[, -1L, drop = FALSE]  # drop intercept
    }
} else {
    prepare_matrix <- function(gene) {
        # clinical only affect the intercept
        X <- cbind(model.matrix(~ clinical), meth[, net[[gene]], drop = FALSE])
        X[, -1L, drop = FALSE]  # drop intercept
    }
}

fit_enet <- function(gene, alpha = seq(0, 1, len = 11L)^2, s = "lambda.1se",
                     family = "gaussian", scale_gex = FALSE) {
    y <- gex[, gene]
    if (family == "gaussian" && isTRUE(scale_gex))
        y <- (y - mean(y)) / sd(y)
    X <- prepare_matrix(gene)
    mod <- tryCatch(cva.glmnet(X, y, alpha = alpha, family = family),
                    error = function(e) e)
    if ("error" %in% class(mod))
        return(mod)

    # retrieve best model
    score <- vapply(mod$modlist, get_cva_score, 0, s = s)
    ix <- which.min(score)
    alpha <- mod$alpha[ix]
    mod <- mod$modlist[[ix]]
    mod$alpha <- alpha
    mod
}

sparse2dt <- function(mod, lambda) {
    mat <- coef(mod, s = lambda)
    dfi <- summary(mat)  # integer rows/cols
    data.table(probe = rownames(mat)[dfi$i], coef = dfi$x)
}

sample_length <- function(x, y) {
    nx <- unique(nchar(x))
    ny <- unique(nchar(y))
    if (length(nx) == 1L && length(ny) == 1L && nx != ny) {
        return(pmin(nx, ny))
    }
    return(NA_integer_)
}


# Read Data --------------------------------------------------------------------

message("Reading data...")

# methylation
meth <- readRDS(meth_file)
if (!all(grepl("^TCGA", rownames(meth))))
    meth <- t(meth)

# gene expression
gex <- readRDS(gex_file)
if (!all(grepl("^TCGA", rownames(meth))))
    gex <- t(gex)

# networks
net <- fread_wrapper(net_file)
if (!is.null(cluster_file)) {
    cluster <- fread_wrapper(cluster_file)
    probe_id_col <- colnames(cluster)[1L]
    if (is.null(cluster_scheme)) {
        cluster_scheme <- colnames(cluster)[2L]
        warning(sprintf("No cluster scheme specified. Using %s by default.", cluster_scheme))
    } else if (!cluster_scheme %in% colnames(cluster)) {
        stop(paste(cluster_scheme, "is not a column of", cluster_file))
    }
    # merge net and clusters and rename
    message("Creating cluster-genes adjacency matrix...")
    net <- net[cluster, on = c("probe" = probe_id_col)][, c("gene", cluster_scheme), with = FALSE]
    net <- unique(net)  # duplicates because multiple probes per cluster
    setnames(net, old = cluster_scheme, new = "probe")
} else {
    setnames(net, 2L, "probe")
}

user_samples <- if (is.null(samples_file)) NULL else scan(samples_file, what = "", sep = "\n")
user_genes   <- if (is.null(genes_file))   NULL else scan(genes_file,   what = "", sep = "\n")


# Filter Data ------------------------------------------------------------------

# intersect samples
message("Intersecting samples across datasets...")
Nsample <- sample_length(rownames(meth), rownames(gex))  # sample character length
if (!is.na(Nsample)) {
    warning(sprintf("Truncating sample names to be %d characters long", Nsample))
    rownames(meth) <- substr(rownames(meth), 1L, Nsample)
    rownames(gex)  <- substr(rownames(gex),  1L, Nsample)
}
samples <- intersect(rownames(meth), rownames(gex))

if (!is.null(user_samples)) {
    samples <- intersect(samples, user_samples)
}

# sanity check
Nsamples <- length(samples)
stopifnot(Nsamples > 2L)
if (Nsamples/nrow(meth) < 0.1 && Nsamples/nrow(gex) < 0.1)
    warning("The common samples of the datasets are less than 10%")

meth <- meth[samples, , drop = FALSE]
gex  <-  gex[samples, , drop = FALSE]
clinical <- factor(substr(samples, 14L, 15L))  # sample_type
cat(samples, file = paste0(prefix, "_samples.txt"), sep = "\n")

if (length(levels(clinical)) == 1L) {
    prepare_matrix <- function(gene) meth[, net[[gene]], drop = FALSE]
}

# user-genes
if (!is.null(user_genes))
    gex <- gex[, intersect(colnames(gex), user_genes), drop = FALSE]

message("Applying filters...")

# variance filters
var_genes <- apply(gex, 2L, mad, na.rm = TRUE)
var_genes <- (var_genes <= gex_min_mad & !is.na(var_genes))
Nvar_genes <- sum(var_genes, na.rm = TRUE)
if (Nvar_genes > 0L)
    warning(paste(Nvar_genes, "genes removed because of low variance (mad)"))
gex <- gex[, !var_genes, drop = FALSE]

# range filters
uq_genes <- apply(gex, 2L, quantile, probs = .75, na.rm = TRUE)
uq_genes <- (uq_genes <= gex_min_uq & !is.na(uq_genes))
Nuq_genes <- sum(uq_genes, na.rm = TRUE)
if (Nuq_genes > 0L)
    warning(paste(Nuq_genes, "genes removed because of low range"))
gex <- gex[, !uq_genes, drop = FALSE]

# convert to integers if family == poisson
if (family == "poisson") {
    # TODO: allow other de-normalization schemes
    warning("De`log`ing count data as (2^x - 1)")
    gex <- round(2^gex) - 1
}

# TODO: impute of NAs
incomp <- apply(is.na(meth) | is.infinite(meth), 2L, any)  # infinities = min/max of NULLs
Nincomp <- sum(incomp)
if (Nincomp > 0L)
    warning(paste(Nincomp, "probes have been removed from the methylation matrix,",
                  "because they contain NAs"))
meth <- meth[, !incomp, drop = FALSE]
meth_clusters <- colnames(meth)
cat(meth_clusters, file = paste0(prefix, "_meth-clusters.txt"), sep = "\n")
colnames(meth) <- paste0("X", 1:ncol(meth))

# clean up network
net <- net[gene %in% colnames(gex) & probe %in% meth_clusters]
net[, probe := factor(probe, meth_clusters, colnames(meth))]
net <- with(net, split(as.character(probe), gene))

# clean-up genes to avoid unnecessary jobs
genes <- intersect(colnames(gex), names(net))
genes <- genes[order(apply(gex[, genes], 2L, mad), decreasing = TRUE)]  # prioritize gene models
gex <- gex[, genes, drop = FALSE]

# Main -------------------------------------------------------------------------


message(sprintf("Training %d models...", length(genes)))
gene_models <- mclapply(genes, fit_enet, s = lambda, family = family, scale_gex = scale_gex)
names(gene_models) <- genes

message("Writing results...")
# write results
saveRDS(gene_models, paste0(prefix, "_enet-models.rds"))  # write errors as well

# deal with failed genes
failed_genes <- vapply(gene_models, function(x) "error" %in% class(x), FALSE)
Nfailed <- sum(failed_genes)
if (Nfailed > 0L) {
    warning(sprintf("%d gene-models failed. See %s_failed-genes.txt", Nfailed, prefix))
    failed_dt <- data.table(
        gene = genes[failed_genes],
        error_msg = vapply(gene_models[failed_genes], "[[", "", "message")
    )
    fwrite(failed_dt, paste0(prefix, "_enet-failed-genes.txt"), sep = "\t")

    gene_models <- gene_models[!failed_genes]
    genes <- names(gene_models)
}

# betas
meth_coefs <- lapply(gene_models, sparse2dt, lambda = lambda)
names(meth_coefs) <- genes
meth_coefs <- rbindlist(meth_coefs, fill = TRUE, idcol = "gene")
meth_coefs[grep("^X", probe), probe := meth_clusters[as.integer(sub("^X", "", probe))]]
if (!is.null(cluster_scheme)) setnames(meth_coefs, "probe", cluster_scheme)
fwrite(meth_coefs, paste0(prefix, "_enet-params.csv"))

# hyperparams
hyper <- cbind(
    data.table(
        gene = genes,
        alpha = vapply(gene_models, "[[", 0, "alpha"),
        lambda = vapply(gene_models, "[[", 0, lambda)
    ),
    t(vapply(gene_models, get_cva_attr, rep(0, 7L)))
)

fwrite(hyper, paste0(prefix, "_enet-hyperparams.csv"))

# plot diagnostics
genes <- hyper[R2 > 0.1][order(-R2)][["gene"]]  # plot for best to worst model ignoring the worst
pdf(file = paste0(prefix, "_enet-diagnostics.pdf"), width = 9, height = 7)
for (gene in genes) {
    mod <- gene_models[[gene]]
    mod$x <- prepare_matrix(gene)
    mod$y <- gex[, gene]
    tryCatch(plotres(mod,
                     info = TRUE, caption = paste(gene, "diagnostics")),
             error = function(e) NULL)
}
dev.off()
