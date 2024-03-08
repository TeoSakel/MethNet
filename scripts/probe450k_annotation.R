#!/usr/bin/env Rscript
doc <- "
Get annotation for Illumina's Human Methylation 450k probes.
Based on hg19 genome.

Usage: probe450k_annotation.R <output>
"


suppressPackageStartupMessages({
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(data.table)
    library(docopt)
})

opts <- docopt(doc)

probes <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probes <- as.data.table(probes, keep.rownames = "probe")
# drop Name (same as rownames, and make implicit FALSE explicit)
probes[, `:=`(Name = NULL,
              DHS = DHS == "TRUE",
              Random_Loci = Random_Loci == "TRUE",
              Methyl27_Loci = Methyl27_Loci == "TRUE",
              Enhancer = Enhancer == "TRUE")]

ext <- tools::file_ext(opts$output)
stopifnot(ext %in% c("csv", "tsv"))
sep <- ifelse(ext == "csv", ",", "\t")
dir.create(dirname(opts$output), showWarnings = FALSE, recursive = TRUE)
fwrite(probes, opts$output, sep=sep)

