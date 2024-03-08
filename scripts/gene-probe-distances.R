#!/usr/bin/env Rscript
doc <- "
Description:
Given annotation files for genes and probes it computes the distance matrices between
gene-tss and probes per chromosome. The results are written as an R-list
object (.rds). If the --tsv flag is given in which case a tsv file is written per chromosome.
Results are written as <prefix>_distances_list.rds or <prefix>_distances_<chr>.tsv.
The first column of <genes> and <probes> is treated as ID.

Usage:
  gene-probe-distances.R [options] <genes> <probes> <prefix>

Options:
  --strand               Ignore strand when computing distances
  --tes                  Compute distance to TES instead of TSS
  --tsv                  Write tab-separated files per chromosome.
  -v --verbose           Write more messages
  --gene-chr=<COL>       Name of chromosome column for gene annotation
  --gene-start=<COL>     Name of start-position column for genes
  --gene-end=<COL>       Name of end-position column for genes
  --probe-chr=<COL>      Name of chromosome column for probes
  --probe-start=<COL>    Name of position column for probes
  --probe-end=<COL>      Name of position column for probes
  --genetype=<type>      Type of genes to filter.
"

# Parameters --------------------------------------------------------------

# hardcoded
chrom_names <- c("seqnames", "seqname", "chromosome", "chrom", "chr", "chromosome_name", "seqid")
start_names <- c("start", "chromStart", "pos")
end_names   <- c("end", "stop", "chromEnd", "chromStop")
chromosomes <- paste0("chr", c(1:22, "X", "Y"))

suppressPackageStartupMessages(library(docopt))
argv <- docopt(doc)

genes_file <- argv[["genes"]]
probes_file <- argv[["probes"]]
prefix <- argv[["prefix"]]

# options
ignore.strand <- argv[["strand"]]
export_tsv    <- argv[["tsv"]]
use_tes       <- argv[["tes"]]
verbose       <- argv[["verbose"]]
gene_chr      <- argv[["gene-chr"]]
gene_start    <- argv[["gene-start"]]
gene_end      <- argv[["gene-end"]]
probe_chr     <- argv[["probe-chr"]]
probe_start   <- argv[["probe-start"]]
probe_end     <- argv[["probe-end"]]
genetype      <- argv[["genetype"]]

outdir <- dirname(prefix)
if (!dir.exists(outdir))
    dir.create(outdir)

# TODO: use a proper logger
say <- ifelse(verbose, message, function(msg) invisible(NULL))


# Load Packages  ----------------------------------------------------------

say("Loading required packages...")
suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(GenomicRanges)
    library(GenomicAlignments)
})

infer_column <- function(column, tab, preset) {
    if (!is.null(column))
        return(column)

    col_matches <- na.omit(match(preset, colnames(tab)))
    colnames(tab)[col_matches[1L]]
}

# Read Coordinates --------------------------------------------------------

## Genes
say("Reading genes' coordinates...")
# genes <- rtracklayer::import(genes_file)
# genes <- subset(genes, seqnames(genes) %in% chromosomes & type == "gene")
# names(genes) <- genes$gene_name
genes <- fread(genes_file)
if (!is.null(genetype))
    genes <- genes[biotype == genetype]
# columns
gene_id <- colnames(genes)[1L]
gene_chr <- infer_column(gene_chr, genes, chrom_names)
gene_start <- infer_column(gene_start, genes, start_names)
gene_end <- infer_column(gene_end, genes, end_names)
# filter
genes <- na.omit(genes[get(gene_chr) %in% chromosomes,
                       c(gene_id, gene_chr, gene_start, gene_end, "strand"),
                       with = FALSE])
# make GRanges
gene_ids <- genes[[1L]]
genes <- makeGRangesFromDataFrame(genes,
                                  keep.extra.columns = FALSE,
                                  seqnames.field = gene_chr,
                                  start.field = gene_start,
                                  end.field = gene_end,
                                  strand.field = "strand")
names(genes) <- gene_ids
genes <- genes[order(genes)]

## Probes
say("Reading probes' coordinates...")
probes <- fread(probes_file)
probe_id <- colnames(probes)[1L]
probe_chr <- infer_column(probe_chr, probes, chrom_names)
probe_start <- infer_column(probe_start, probes, start_names)
probe_end <- infer_column(probe_end, probes, end_names)
if (is.na(probe_end)) {
    probes[, end := get(probe_start) + fifelse(strand == "+", 1L, -1L)]
    probe_end <- "end"
}

# filter
probes <- na.omit(probes[get(probe_chr) %in% chromosomes,
                         c(probe_id, probe_chr, probe_start, probe_end, "strand"),
                         with = FALSE])
# make GRanges
probe_ids <- probes[[1L]]
probes <- makeGRangesFromDataFrame(as.data.frame(probes),
                                   keep.extra.columns = FALSE,
                                   seqnames.field = probe_chr,
                                   start.field = probe_start,
                                   end.field = probe_end,
                                   ignore.strand = FALSE)
names(probes) <- probe_ids
probes <- probes[order(probes)]


# Compute Distances -------------------------------------------------------

# Overlaps
say("Finding overlaps...")
overlaps <- findOverlaps(genes, probes, ignore.strand = TRUE)
overlaps <- data.table(gene = names(genes)[from(overlaps)],
                       probe = names(probes)[to(overlaps)])

# Distances
compute_dist <- function(chr, ignore.strand = FALSE, x0 = "tss") {
    chr_genes <- genes[seqnames(genes) == chr]
    if (x0 == "tss") {
        x <- start(resize(chr_genes, 1L, fix = "start"))
    } else if (x0 == "tes") {
        x <- start(resize(chr_genes, 1L, fix = "end"))
    }
    chr_probes <- probes[seqnames(probes) == chr]
    prb <- start(chr_probes)
    d <- -1L * outer(x, prb, FUN = '-')  # positive == downstream (5-3 dir)
    rownames(d) <- names(chr_genes)
    colnames(d) <- names(chr_probes)
    if (!ignore.strand) {
        mul <- ifelse(strand(chr_genes) == "-", -1L, 1L)  # positive == parallel to transcription
        d <- d * mul  # multiply by row
    }
    d
}

say("Computing distances...")
x0 <- ifelse(use_tes, "tes", "tss")
dmat <- lapply(chromosomes, compute_dist, ignore.strand = ignore.strand, x0 = x0)
names(dmat) <- chromosomes

dmat[["overlaps"]] <- overlaps


# Write results -----------------------------------------------------------

say("Writing list...")
saveRDS(dmat, paste0(prefix, "_distances_list.rds"))

if (export_tsv) {
    say("Writing tsv-files...")
    for (chr in chromosomes) {
        fwrite(dmat[[chr]], paste0(prefix, "_distances_", chr, ".tsv"), sep = "\t")
    }
    fwrite(overlaps, paste0(prefix, "_overlaps.tsv"), sep = "\t")
}

say("Done!")

