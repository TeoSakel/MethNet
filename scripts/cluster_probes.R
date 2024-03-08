#!/usr/bin/env Rscript
doc <- "
Clusters CpG probes based on their DNA coordinate distance and a given window size.

It produces 2 types of clustering:
    - single: no probes belonging to different clusters have distance greater than the window size
    - complete: no probes within the same cluster have distance greater than the window size

Usage: cluster_probes.R <probe_coords> <window_size> <output>

Arguments:
    probe_coords    A comma/tab separated file with:
                        1. a `probe_id` as it's 1st column
                        2. a `chromosome` column named as chrom/chr/etc
                        3. a `position` column name as pos/start/chromStart
    window_size     The cluster radius/diameter in bp
    output          Path to write results
"

library(docopt)
argv <- docopt(doc)
# argv <- docopt(doc, c("../data/annot/methyl450k_hg19_pancan.tsv", "200", "test200.csv"))

probe_file <- argv$probe_coords
output_file <- argv$output
window_size <- suppressWarnings(as.integer(argv$window_size))

# Checks
if (!file.exists(probe_file))
    stop(sprintf("Input file '%s' does not exists", probe_file))

if (is.na(window_size))
    stop(sprintf("Failed to transform %s to integer", argv$window_size))

chrom_names <- c("seqnames", "seqname", "chromosome", "chrom", "chr", "chromosome_name", "seqid")
start_names <- c("pos", "start", "chromStart")
chromosomes <- paste0("chr", c(1:22, "X", "Y"))


# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

complete_clusters <- function(pos, window_size) {
    # returns a complete clustering of pos (max distance within cluster = window_size)
    # assuming pos is ordered vector of coordinates
    n <- length(pos)
    if (n < 3L || pos[n] - pos[1L] <= window_size)
        return(rep(1L, n))

    outer(pos, pos, FUN = "-") %>%   # very wasteful step (only lower-tri required)
        abs() %>%
        as.dist() %>%
        hclust(method = "complete") %>%
        cutree(h = window_size)
}

if (str_detect(probe_file, "(tsv|tsv.gz)$")) {
    read_probes <- read_tsv
} else if (str_detect(probe_file, "(csv|csv.gz)$")) {
    read_probes <- read_csv
} else {
    stop(paste("Unknown file extention for", probe_file))
}

if (str_detect(output_file, "(tsv|tsv.gz)$")) {
    write_result <- write_tsv
} else if (str_detect(output_file, "(csv|csv.gz)$")) {
    write_result <- write_csv
} else {
    stop(paste("Unknown file extention for", output_file))
}

cluster_range <- function(pos) paste0(pos[1L], "-", pos[length(pos)] + 1L)

infer_column <- function(preset, colNames) {
    col_matches <- na.omit(match(preset, colNames))
    colNames[col_matches[1L]]
}

# Main --------------------------------------------------------------------

window_bp = paste0("_", window_size, "bp")  # to hardcode this info in the results

probes <- suppressMessages(read_probes(probe_file, na = c("", "*", ".")))
probe_id <- sym(colnames(probes)[1L])
chrom <- sym(infer_column(chrom_names, colnames(probes)))
pos <- sym(infer_column(start_names, colnames(probes)))

probes <- probes %>%
    select(!!!c(probe_id, chrom, pos)) %>%
    drop_na() %>%
    filter(str_detect(!!probe_id, "^cg")) %>%  # ignore ^rs probes
    arrange(!!chrom, !!pos) %>%  # crucial assumption
    group_by(!!chrom) %>%
    mutate(adj_dist = c(window_size + 1L, diff(!!pos)),  # window_size + 1L = force new cluster
           single = cumsum(adj_dist > window_size)) %>%       # when true a new cluster is generated
    group_by(!!chrom, single) %>%
    # only need to consider sub-clusters (the others are definetely > window_size)
    mutate(single_range = cluster_range(!!pos),
           complete = complete_clusters(!!pos, window_size)) %>%
    # group_by(!!chrom) %>%
    ungroup() %>%
    mutate(complete = paste(single, complete, sep = ".")) %>%  # make names unique
    group_by(!!chrom, complete) %>%
    mutate(complete_range = cluster_range(!!pos)) %>%
    ungroup() %>%
    transmute(probe = !!probe_id, # to survive transmute
              !!paste0("single", window_bp)   := paste(!!chrom, single_range, sep = ":"),
              !!paste0("complete", window_bp) := paste(!!chrom, complete_range, sep = ":"))

write_result(probes, output_file)

