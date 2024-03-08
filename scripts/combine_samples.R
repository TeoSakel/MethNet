#!/usr/bin/env Rscript
# Parse inputs -----------------------------------------------------------------
doc <- "
Concatanate martrices by rows

Usage: combine_samples.R <output> <matrix>...
"

library(docopt)
argv <- docopt(doc)

mat_files <- argv[["matrix"]]
output <- argv[["output"]]

# Main  ----------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

read_mat <- function(path) {
    mat <- if (tools::file_ext(path) == "rds") {
        mat <- readRDS(path)
        if (!all(grepl("^TCGA", rownames(mat))))
            mat <- t(mat)
        mat <- as.data.table(mat, keep.rownames = "sampleID")
    } else {  # flat file
        mat <- fread(path)
        if (! "sampleID" %in% names(mat))
            mat <- transpose(mat, keep.names = "sampleID", make.names = 1L)
    }
    mat
}

mat <- rbindlist(lapply(mat_files, read_mat), fill = TRUE, use.names = TRUE)

ext <- tools::file_ext(output)
if (ext == "rds") {
    mat <- as.matrix(mat, rownames = "sampleID")
    saveRDS(mat, output)
} else {
    sep <- ifelse(ext == "tsv", "\t", ",")
    fwrite(mat, output, sep = sep)
}

