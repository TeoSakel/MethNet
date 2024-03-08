#!/usr/bin/env Rscript
doc <- "
Augment xenahub gene-annotation (hg19) using Ensembl.

Usage: augment_genes.R <input> <output>
"

library(docopt)
opt <- docopt(doc)

suppressPackageStartupMessages({
    library(EnsDb.Hsapiens.v75)  # GRCh37 == hg19
    library(data.table)
})

genes <- fread(opt$input, select=1:6)
genes[, gene := NULL]
genes <- na.omit(genes)
genes[, biotype := mapIds(EnsDb.Hsapiens.v75, keys = id, keytype = "GENENAME", column = "GENEBIOTYPE")]
fwrite(genes, opt$output, sep = "\t")
