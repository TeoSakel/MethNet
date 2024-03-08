#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(broom)
    library(survival)
    library(data.table)
    library(parallel)
})

Ncores <- future::availableCores()
setDTthreads(Ncores)

meth_files <- snakemake@input[["meth"]]
clin_file <- snakemake@input[["survival"]]
output_coefs <- snakemake@output[["coefs"]]
output_stats <- snakemake@output[["stats"]]
cre <- snakemake@input[["cre"]]
resdir <- dirname(output_coefs)

if (!dir.exists(resdir)) dir.create(resdir, recursive = TRUE)

read_meth <- function(path, cre = NULL) {
    ext <- tools::file_ext(path)
    if (ext == "rds") {
        meth <- readRDS(path)
        if (!is.null(cre))
            meth <- meth[, intersect(cre, colnames(meth)), drop = FALSE]
        meth <- as.data.table(meth, keep.rownames = "sampleID")
        setkey(meth, "sampleID")
    } else if (ext %in% c("csv", "tsv")) {
        if (!is.null(cre)) cre <- c("sampleID", cre)
        meth <- fread(path, select = cre, key = "sampleID")
    } else {
        stop(sprintf("Uknown file extension (%s) for path: %s", ext, path))
    }
    meth
}

fit_survival <- function(x, y, cancer) {
    few_samples <- names(which(table(cancer[!is.na(x)]) <= 10))
    x[cancer %in% few_samples] <- NA_real_
    ix <- which(!is.na(x))
    if (length(ix) < 100) return(NULL)
    fml <- if (length(unique(cancer[ix])) > 1) y ~ x*strata(cancer) else y ~ x
    tryCatch(coxph(fml), error = function(e) {message(e); return(NULL)})
}


# Main ------------------------------------------------------------------------

if (!is.null(cre)) cre <- scan(cre, what = character(), quiet = TRUE)

meth <- lapply(meth_files, read_meth, cre = cre)
meth <- rbindlist(meth, use.names = TRUE, fill = TRUE)
setkey(meth, sampleID)
if (is.null(cre)) cre <- names(meth)[-1]

## Read & cleanup clinical
clin_cols <- c("sample", "cancer type abbreviation", "OS", "OS.time")
clin <- fread(clin_file, select = clin_cols)
clin <- na.omit(clin)
setnames(clin, c("cancer type abbreviation", "sample"), c("cancer", "sampleID"))
clin <- clin[grep("0[1-7]$", sampleID)][
             , submitter_id := substr(sampleID, 1L, 12L) ][
             , if (.N == 1L) .SD, submitter_id ][
             grep("0(1|3|6)", sampleID) ][
             , submitter_id := NULL]
clin <- clin[, if (.N > 10) .SD else NULL, cancer]
setkey(clin, sampleID)

sample_ids <- sort(intersect(clin$sampleID, meth$sampleID))
meth <- meth[sample_ids]
clin <- clin[sample_ids]
meth[clin, `:=`(OS = i.OS, OS.time = i.OS.time, cancer = i.cancer)]
names(cre) <- paste0("X", 1:length(cre))
setnames(meth, cre, names(cre))

y <- meth[, Surv(OS.time, OS)]
setDTthreads(1L)
meth <- mclapply(meth[, names(cre), with = FALSE], fit_survival,
                 y = y, cancer = meth[["cancer"]], mc.cores = Ncores)
names(meth) <- cre
meth <- meth[vapply(meth, class, "") == "coxph"]

setDTthreads(Ncores)
params <- rbindlist(lapply(meth, tidy), idcol="CRE")
params[, term := sub("^x", "CpG", term)]
fwrite(params, output_coefs)

qual <- rbindlist(lapply(meth, glance), idcol="CRE")
fwrite(qual, output_stats)

