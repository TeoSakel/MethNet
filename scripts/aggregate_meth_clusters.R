#!/usr/bin/env Rscript
# Parse inputs -----------------------------------------------------------------
doc <- "
Aggregates methylation values according to clusters.
Aggregation methods available are: mean, max, uq (upper quantile)

Arguments:
    - meth_matrix: path to methylation matrix (rows = probes, cols = samples)
    - clusters: path to cluster matrix (first cols must be `probe_id` the rest clustering-schemes)
    - scheme: which column of the <clusters> file to use for clustering
    - output: path to write results

Usage:
    aggregate_meth_clusters.R [options] mean <meth_matrix> <clusters> <scheme> <output>
    aggregate_meth_clusters.R [options] max  <meth_matrix> <clusters> <scheme> <output>
    aggregate_meth_clusters.R [options] uq   <meth_matrix> <clusters> <scheme> <output>

Options:
    -z, --compress       Compress resulting file (if rds)
    -p ID --probes=ID    probe id column in <clusters> [default: probe]
    -L LN --length=LN    Length to truncate samples if not positive it skips truncation [default: 15]
"

library(docopt)

argv <- docopt(doc)

meth_file      <- argv[["meth_matrix"]]
cluster_file   <- argv[["clusters"]]
cluster_scheme <- argv[["scheme"]]
output_file    <- argv[["output"]]
probe_id       <- argv[["probes"]]
L              <- as.integer(argv[["length"]])
compress       <- argv[["compress"]]

stopifnot(grepl("(tsv|csv|rds)$", output_file))  # crash early
stopifnot(is.integer(L))


# Library ----------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
setDTthreads(percent=100)

agg.fun <- if (argv$mean) {
    mean
} else if (argv$max) {
    max
} else if (argv$uq) {
    function(x, ...) ifelse(length(x) < 4L, max(x, ...), unname(quantile(x, probs = 0.75, ...)))
} else {
    stop("Unknown aggregation function")
}


# Main -------------------------------------------------------------------------

# clusters
clusters <- fread(cluster_file, select = c(probe_id, cluster_scheme), key = probe_id)

# methylation matrix
if (tools::file_ext(meth_file) == "rds") {  # already a matrix
    meth <- as.data.table(readRDS(meth_file), keep.rownames = probe_id)
} else {  # flat file
    meth <- fread(meth_file)
    setnames(meth, old = 1L, new = probe_id)
}
meth <- na.omit(meth)
setkeyv(meth, probe_id)

# map to clusters and aggregate
samples <- colnames(meth)[-1L]
meth <- clusters[meth, nomatch = 0L ][  # inner join
    , if (.N > 1) lapply(.SD, agg.fun, na.rm = TRUE) else .SD, by = cluster_scheme, .SDcols = samples]
setkeyv(meth, cluster_scheme)

# merge samples
meth <- transpose(meth, keep.names = "sampleID", make.names = cluster_scheme)
L <- ifelse(L > 0, L, meth[, max(nchar(sampleID))])
meth[, sampleID := substr(sampleID, 1L, L)]
meth <- meth[, if (.N > 1) lapply(.SD, mean, na.rm = TRUE) else .SD, sampleID]

# write results
file_ext <- tools::file_ext(output_file)
if (file_ext == "rds") {
    # R-matrix columns = cpg clusters, rows = samples
    meth <- as.matrix(meth, rownames = "sampleID")
    saveRDS(meth, output_file, compress = compress)
} else {
    fwrite(meth, output_file, sep = ifelse(file_ext == "csv", ",", "\t"))
}

