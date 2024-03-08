#!/usr/bin/env Rscript
doc <- "
Description:
Given an R-list (rds file) with chr as keys and distance between-matrices
it generates an edge-list with columns: chr, gene, probe, and distance.
In addition if 'overlaps' exists as a key or suffix to a file,
an extra column 'overlap' is added.

Usage: gene-probe-network.R [options] <distances> <genes> <output> <cutoff>

Options:
  -c, --combine    Combine results
  --ref=REF        Reference point [default: tss]
"
suppressPackageStartupMessages(library(docopt))
argv <- docopt(doc)

dist_file <- argv[["distances"]]
cutoff <- suppressWarnings(as.numeric(argv[["cutoff"]]))
output <- argv[["output"]]
gene_file <- argv[["genes"]]
refpoint <- tolower(argv[["ref"]])
stopifnot(refpoint %in% c("tss", "body"))

if (is.na(cutoff))
    stop(paste("Cutoff must be numeric.", argv[["cutoff"]], "was given"))

dmat <- readRDS(dist_file)
chromosomes <- setdiff(names(dmat), "overlaps")


# Libraries  ---------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
genes <- fread(gene_file, select=1:6)
setnames(genes, c("id", "chrom", "start", "end", "strand", "biotype"))
gene_width <- genes[, setNames(abs(end - start), id)]
overlaps <- as.data.table(dmat[["overlaps"]])

mat2dt <- function(mat, cutoff, ref = "tss") {
    ig <- intersect(rownames(mat), names(gene_width))
    mat <- mat[ig, , drop = FALSE]
    if (ref == "tss") {
        gw <- 0L
    } else if (ref == "body") {
        gw <- gene_width[ig]
    } else {
        stop(sprintf("Uknown reference point %s", ref))
    }
    Nr <- nrow(mat)
    indx <- as.data.table(which(mat >= -cutoff & mat <= gw + cutoff, arr.ind = TRUE))  # subscript index
    if (nrow(indx) == 0)
        return(NULL)
    indx[, ':='(gene  = rownames(mat)[row],
                probe = colnames(mat)[col],
                dist  = mat[(col - 1L) * Nr + row])]  # linear index
    indx[, dist_tes := dist - gene_width[gene]]
    return(indx[ , .(gene, probe, dist, dist_tes)])
}

merge_overlaps <- function(indx, overlaps) {
    if (!nrow(overlaps))
        return(indx)
    indx[overlaps, on = c("gene", "probe"), overlap := TRUE]  # assign TRUE to matching rows
    indx[, overlap := !is.na(overlap)]
    return(indx)
}

# Main --------------------------------------------------------------------

if (argv[["combine"]]) {
    ddt <- rbindlist(lapply(dmat[chromosomes], mat2dt, cutoff = cutoff, ref = refpoint))
    ddt <- merge_overlaps(ddt, overlaps)
    ext <- tools::file_ext(output)
    if (ext != "tsv") output <- paste0(output, ".tsv")
    fwrite(ddt, output, sep = "\t")
} else {
    for (chr in chromosomes) {
        ddt <- mat2dt(dmat[[chr]], cutoff, ref = refpoint)
        ddt <- merge_overlaps(ddt, overlaps)
        fwrite(ddt, paste0(output, "_", chr, ".tsv"), sep = "\t")
    }
}
