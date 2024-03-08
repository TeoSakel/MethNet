#!/usr/bin/env Rscript
doc <- "
Compute summary statistic of methylation status of probes.

Usage: meth_cluster_beta_summary.R <output> <meth>...
"

library(docopt)
argv <- docopt(doc)
output <- argv[["output"]]
meth_files <- argv[["meth"]]
names(meth_files) <- tools::file_path_sans_ext(basename(meth_files))
Ncores <- future::availableCores()


suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(parallel)
})


summarize_mat <- function(mat)  data.table(cluster_id = rownames(mat),
                                           beta_mean = rowMeans2(mat, na.rm = T),
                                           beta_sd = rowSds(mat, na.rm = T),
                                           Nsamples = rowSums(!is.na(mat)),
                                           key = "cluster_id")

read_mat2dt <- function(path) {
    mat <- readRDS(path)
    if (!all(grepl("^TCGA", colnames(mat)))) mat <- t(mat)
    normal <- grep("-11[A-Z]?$", colnames(mat))
    if (length(normal) > 0) {
        df_tumor <- summarize_mat(mat[, -normal, drop = FALSE])
        df_normal <- summarize_mat(mat[, normal, drop = FALSE])
        df <- merge(df_tumor, df_normal, key = "cluster_id", suffixes = c("_tumor", "_normal"))
    } else {
        df <- summarize_mat(mat)
        setnames(df, 2:4, function(x) paste0(x, "_tumor"))
    }
    df
}

setDTthreads(1L)
DT <- mclapply(meth_files, read_mat2dt, mc.cores = Ncores)
setDTthreads(Ncores)
DT <- rbindlist(DT, use.names = TRUE, fill = TRUE, idcol = "cancer")
fwrite(DT, output)
