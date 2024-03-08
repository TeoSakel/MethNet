#!/usr/bin/env Rscript
# Parse Inputs  ----------------------------------------------------------------
doc <- "
Find CpG-clusters overlaping genes

Usage: gene_cluster_overlap.R  <probe_network> <cluster_file> <cluster_scheme> <outfile>
"

suppressPackageStartupMessages({
    library(docopt)
    library(data.table)
})

argv <- docopt(doc)

get_overlaps <- function(net_file, cluster_file, cluster_scheme) {
    net <- fread(net_file)
    if (!cluster_scheme %in% colnames(net)) {
        clusters <- fread(cluster_file, select = c("probe", cluster_scheme))
        net <- net[clusters, nomatch = 0, on = "probe"]  # merge net + clusters -> map clusters to genes
    }

    net[, .(overlap = any(overlap)), c("gene", cluster_scheme)]
}

net <- with(argv, get_overlaps(probe_network, cluster_file, cluster_scheme))
fwrite(net, argv[["outfile"]])
