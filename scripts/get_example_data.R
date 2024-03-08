#!/usr/bin/env Rscript
# Parse Inputs  ----------------------------------------------------------------
doc <- "
Combine training data to visualize relationships.

Usage: get_example_data.R <meth_dir> <rna> <output> (<gene> <elem>)...
"

library(docopt)
argv <- docopt(doc)

meth_dir <- argv[["meth_dir"]]
rna_file <- argv[["rna"]]
genes <- argv[["gene"]]
elems <- argv[["elem"]]
output <- argv[["output"]]


# Libraries  -------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

read_and_subset <- function(path, elem) {
    mat <- readRDS(path)
    if (all(elem %in% rownames(mat))) {
        mat <- t(mat[elem, , drop = FALSE])
    } else if (all(elem %in% colnames(mat))) {
        mat <- mat[, elem, drop = FALSE]
    } else {
        xrow <- setdiff(elem, rownames(mat))
        xcol <- setdiff(elem, colnames(mat))
        if (length(xrow) < length(xcol)) {
            stop(sprintf("Some elements not in matrix %s rows: [%s]",
                         path, paste(xrow, collapse = ", ")))
        } else {
            stop(sprintf("Some elements not in matrix %s colums: [%s]",
                         path, paste(xcol, collapse = ", ")))
        }
    }
    as.data.table(mat, keep.rownames = "sampleID")
}


# Main  -------------------------------------------------------------------------

setDTthreads(percent=100)

associations <- data.table(elem = elems, gene = genes)

meth_files <- list.files(meth_dir, full.names = TRUE)
meth <- lapply(meth_files, read_and_subset, elem = elems)
names(meth) <- sub(".rds", "", basename(meth_files))
meth <- rbindlist(meth, idcol = "cancer")
meth[, sampleID := substr(sampleID, 1L, 15L)]
meth <- meth[, lapply(.SD, mean, na.rm = TRUE), .(cancer, sampleID)]
meth <- melt(meth, id.vars = c("cancer", "sampleID"), variable.name = "elem", value.name = "beta")
meth[associations, gene := i.gene, on = "elem"]

rna <- read_and_subset(rna_file, genes)
rna <- melt(rna, id.vars = "sampleID", variable.name = "gene", value.name = "expr")

DT <- merge(meth, rna, by = c("sampleID", "gene"), all = FALSE)
sep <- ifelse(tools::file_ext(output) == "tsv", "\t", ",")
fwrite(DT, output, sep = sep)

