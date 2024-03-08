#!/usr/bin/env Rscript
# Parse inputs -----------------------------------------------------------------

doc <- "
Split pancancer matrix to individual cancer matrices.

Arguments:
    - matrix  : path to pancancer matrix (rows = probes, cols = samples)
    - clinical: path to clinical data files
    - outdir  : path to write results

Usage:
    split_pancancer.R [--compress] <matrix> <clinical> <outdir>
"

library(docopt)

argv <- docopt(doc)

compress <- argv[["compress"]]
mat      <- argv[["matrix"]]
clinical <- argv[["clinical"]]
outdir   <- argv[["outdir"]]

rds <- sub(".(tsv|csv)$", ".rds", mat)


# Main -------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

meta <- fread(clinical, select = c(1, 3), col.names = c("sampleID", "cancer"), key="sampleID")
mat <- fread(mat)

mat <- as.matrix(mat, rownames = 1L)
saveRDS(mat, rds, compress = compress)

mat <- t(mat)  # samples x genes
mat <- as.data.table(mat, keep.rownames = "sampleID")
setkey(mat, "sampleID")

mat[meta, cancer := tolower(i.cancer)]
mat <- split(mat[!is.na(cancer)], by = "cancer", keep.by = FALSE)
mat <- lapply(mat, as.matrix, rownames = 1L)
for (cancer in names(mat))
    saveRDS(mat[[cancer]], file.path(outdir, paste0(cancer, ".rds")), compress = compress)

