#!/usr/bin/env Rscript
# Parse Inputs ------------------------------------------------------------

doc <-"
Usage: gene-cluster-network.R [options] <net> <clusters> <scheme> <output>

Options:
  --genes=<path>    path to gene annotations to compute relative distance
"

library(docopt)
argv <- docopt(doc)

net_file <- argv[["net"]]
cluster_file <- argv[["clusters"]]
cluster_scheme <- argv[["scheme"]]
output <- argv[["output"]]
gene_file <- argv[["genes"]]


# Main --------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
setDTthreads(percent=100)

# read data
net <- fread(net_file)
clusters <- fread(cluster_file)

# check cluster scheme
probe_id_col <- colnames(clusters)[1L]
if (is.null(cluster_scheme)) {
    cluster_scheme <- colnames(clusters)[2L]
    warning(sprintf("No cluster scheme specified. Using %s by default.", cluster_scheme))
} else if (!cluster_scheme %in% colnames(clusters)) {
    stop(paste(cluster_scheme, "is not a column of", cluster_file))
}

# merge net and clusters and rename
message("Creating cluster-genes adjacency matrix...")
clusters <- clusters[, c(probe_id_col, cluster_scheme), with = FALSE]
setnames(clusters, c("probe", "cluster"))
net[clusters, cluster := i.cluster, on = "probe"]
net <- net[!is.na(cluster),
           .(probes = .N, overlap = sum(overlap),
             dist_min     = min(dist),     dist_avg     = mean(dist),     dist_max     = max(dist),
             dist_min_tes = min(dist_tes), dist_avg_tes = mean(dist_tes), dist_max_tes = max(dist_tes)),
           .(gene, cluster)]

# compute relative distance
if (!is.null(gene_file)) {
    genes <- fread(gene_file, select = c("id", "chromStart", "chromEnd"))
    genes[, width := abs(chromEnd - chromStart)]
    net[genes, width := i.width, on = c("gene" = "id")]
    net[, `:=`(reldist_min = dist_min / width,
               reldist_max = dist_max / width,
               reldist_avg = dist_avg / width,
               width = NULL)]
}

# write results
dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
sep <- ifelse(tools::file_ext(output) == "tsv", "\t", ",")
fwrite(net, output, sep = sep)
