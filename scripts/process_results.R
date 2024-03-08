#!/usr/bin/env Rscript
# Parse Input -------------------------------------------------------------
doc <- "

Usage: process_results.R [options] <net> <hyper> <coef> <samples> <outdir>

Options:
  --min_samples=N      minimum number of cluster per samples [default: 100]
  --min_coef=VAL       minimum absolute value of coefficient [default: 0.0001]
  --promoter_up=BP     promoter border upstream [default: 2000]
  --promoter_down=BP   promoter border upstream [default: 200]
"

library(docopt)

argv <- docopt(doc)

net_file    <- argv[["net"]]
hyper_file  <- argv[["hyper"]]
coef_file   <- argv[["coef"]]
sample_file <- argv[["samples"]]
outdir      <- argv[["outdir"]]
min_samples <- as.integer(argv[["min_samples"]])
min_coef    <- as.numeric(argv[["min_coef"]])
prom_up     <- as.integer(argv[["promoter_up"]])
prom_dn     <- as.integer(argv[["promoter_down"]])


# Main --------------------------------------------------------------------

suppressPackageStartupMessages({
    library(data.table)
    library(mgcv)
})
setDTthreads(percent=100)


# Read Data
sample_cancer <- fread(sample_file)
sample_counts <- sample_cancer[, .(Nsamples = .N), cancer][Nsamples >= min_samples]

net <- fread(net_file)

hyper <- fread(hyper_file)
hyper <- sample_counts[hyper, on = "cancer", nomatch = NULL]
hyper[, nvar := nzero + nnzero]

params <- fread(coef_file)
params <- params[cancer %in% sample_counts$cancer & abs(coef) > min_coef]
params[net, dist := i.dist_avg, on = c("gene", "probe" = "cluster")]

# All Associations
associations <- params[!is.na(dist)]
setnames(associations, "probe", "cluster")
params <- params[is.na(dist)]
setnames(params, "probe", "variable")
params[, dist := NULL]
fwrite(params, file.path(outdir, "intercepts.csv"))

associations[, contrib := abs(coef) / sum(abs(coef)), c("cancer", "gene")]
associations[hyper, score := log2(contrib) + log2(i.nvar), on = c("cancer", "gene")]
associations[, dist := NULL]
fwrite(associations, file.path(outdir, "associations.csv"))

# MethNet
hyper[, W := R2/sum(R2), gene]
associations[hyper, w_coef := coef * i.W, on = c("cancer", "gene")]
methnet <- associations[grep("chrY", cluster, invert = TRUE),
                        .(coef = sum(w_coef), nup = sum(coef > 0), ndn = sum(coef < 0)),
                        by = c("gene", "cluster")]
methnet[net, c("dist", "overlap") := .(i.dist_avg, i.overlap), on = c("gene", "cluster")]
setcolorder(methnet, c("gene", "cluster", "dist", "overlap"))
methnet[hyper[, .N, gene], ncancer := i.N, on = "gene"]
methnet[, contrib := abs(coef) / sum(abs(coef)), by = "gene"]
methnet[net[, .N, gene], score := log2(contrib) + log2(i.N), on = "gene"]
methnet[, promoter := dist >= -prom_dn & dist <= prom_up]
methnet[, intergenic := !promoter & overlap == 0]
fwrite(methnet, file.path(outdir, "methnet.csv"))

# Cluster Scores
cluster_score <- methnet[intergenic == TRUE]
cluster_score[, score := resid(gam(score ~ s(abs(dist), bs = "ts"), weights = ncancer), type = "pearson")]
cluster_score <- cluster_score[, .(score = sum(score), igain = sum(pmax(score, 0)),
                                   ninter = .N, nup = sum(coef > 0), ndn = sum(coef < 0)),
                               cluster]
cluster_score[, score_rank := rank(score) / .N]
fwrite(cluster_score, file.path(outdir, "cluster_score.csv"))
