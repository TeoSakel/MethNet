#!/usr/bin/env Rscript
doc <- "
Combine results of fit_glmnet_models.R across multiple cancers.

Usage: combine_enet_models.R <model_dir> <prefix>
"

library(docopt)
opts <- docopt(doc)

suppressPackageStartupMessages({
    library(tidyverse)
})

model_dir <- opts[["model_dir"]]
prefix <- opts[["prefix"]]
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)

## Params
params_suffix =  "_enet-params.csv"
outfile <- sprintf("%s_params.csv", prefix)
params <- list.files(model_dir, pattern = params_suffix, recursive=FALSE)
file.path(model_dir, params) %>%
    set_names(str_remove(params, params_suffix)) %>%
    map_dfr(read_csv, col_types = "ccd", .id = "cancer") %>%
    write_csv(outfile)

## Hyperparams
hyper_suffix <- "_enet-hyperparams.csv"
outfile <- sprintf("%s_hyperparams.csv", prefix)
hyper <- list.files(model_dir, pattern = hyper_suffix, recursive=FALSE)
file.path(model_dir, hyper) %>%
    set_names(str_remove(hyper, hyper_suffix)) %>%
    map_dfr(read_csv, col_types = "cdddddddii", .id = "cancer") %>%
    write_csv(outfile)

# Meth cluster
cluster_suffix <- "_meth-clusters.txt"
outfile <- sprintf("%s_clusters.csv", prefix)
cluster <- list.files(model_dir, pattern = cluster_suffix)
cluster_tbl <- map(file.path(model_dir, cluster), scan, what = "")
tibble(cancer = rep(str_remove(cluster, cluster_suffix), map_int(cluster_tbl, length)),
       cluster = unlist(cluster_tbl)) %>%
    write_csv(outfile)

# Samples
sample_suffix <- "_samples.txt"
outfile <- sprintf("%s_samples.csv", prefix)
sample_files <- list.files(model_dir, pattern = sample_suffix)
sample_tbl <- map(file.path(model_dir, sample_files), scan, what = "")
tibble(cancer = rep(str_remove(sample_files, sample_suffix), map_int(sample_tbl, length)),
       submitter_id = unlist(sample_tbl)) %>%
    write_csv(outfile)
