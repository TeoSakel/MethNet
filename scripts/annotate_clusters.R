#!/usr/bin/env Rscript
doc <- "
Download genome annotation and overlap with CpG clusters.

Usage: annotate_clusters.R <probe_annot> <cluster_file> <cluster_scheme> <outdir>
"
library(docopt)

argv <- docopt(doc)
probe_annot_file <- argv$probe_annot
cluster_file <- argv$cluster_file
cluster_scheme <- argv$cluster_scheme
outdir <- argv$outdir

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicAlignments)
    library(InteractionSet)
    library(AnnotationHub)
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# Get Cluster Coordinates -------------------------------------------------

probe2cluster <- fread(cluster_file, select = c("probe", cluster_scheme))
setnames(probe2cluster, c("probe", "cluster"))
clusters <- probe2cluster[["cluster"]] |> unique() |> GRanges()


# Aggregate probe annotation  ---------------------------------------------

probe_annot <- fread(probe_annot_file)
probe_annot[probe2cluster, cluster := i.cluster, on = "probe"]
cluster_annot <- probe_annot[
    !is.na(cluster),
    .(Nprobes = .N,
      SNPs = sum(!is.na(Probe_rs)),
      Relation_to_Island = fcase(any(grepl("Shelf", Relation_to_Island)), "Shelf",
                                 any(grepl("Shore", Relation_to_Island)), "Shore",
                                 any(Relation_to_Island == "Island"), "Island",
                                 default = "OpenSea"),
      Enhancer = any(Enhancer),
      DHS = any(DHS),
      Regulatory_Feature = any(!is.na(Regulatory_Feature_Group) & !grepl("Unclassified", Regulatory_Feature_Group)),
      UCSC_RefGene = any(!is.na(UCSC_RefGene_Name))),
    cluster
]
setkey(cluster_annot, cluster)
# fwrite(cluster_annot, file.path(outdir, "meth-cluster_manifest.csv"))

# Encode Segmentation -----------------------------------------------------

cells <- c("Gm12878", "H1hesc", "Hepg2", "Hmec", "Hsmm",
           "Huvec", "K562", "Nhek", "Nhlf")

states <- c("1_Active_Promoter", "2_Weak_Promoter", "3_Poised_Promoter",
            "4_Strong_Enhancer", "5_Strong_Enhancer", "6_Weak_Enhancer",
            "7_Weak_Enhancer", "8_Insulator", "9_Txn_Transition", "10_Txn_Elongation",
            "11_Weak_Txn", "12_Repressed", "13_Heterochrom/lo", "14_Repetitive/CNV",
            "15_Repetitive/CNV")

meta_state_map <- c("Promoter", "Promoter", "Poised_Promoter",
                    "Enhancer", "Enhancer", "Weak_Enhancer", "Weak_Enhancer",
                    "Insulator",
                    "Txn", "Txn", "Txn",
                    "Repressed",
                    "Low_Signal", "Low_Signal", "Low_Signal")

meta_states <- c("Promoter", "Poised_Promoter", "Enhancer", "Weak_Enhancer",
                 "Insulator", "Repressed", "Txn", "Low_Signal")

base_dir <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm"
hmm <- file.path(base_dir, sprintf("wgEncodeBroadHmm%sHMM.bed.gz", cells))
hmm <- lapply(hmm, rtracklayer::import)
names(hmm) <- cells
cells <- Rle(cells, sapply(hmm, length))
hmm <- Reduce(c, hmm)
hmm$cell <- cells
hmm$name <- factor(hmm$name, states)
hmm$score <- NULL
hmm$itemRgb <- NULL
hmm$thick <- NULL
hmm <- sort(hmm)
hmm$state <- factor(meta_state_map[as.integer(hmm$name)], meta_states)

gr <- disjoin(hmm)
hits <- findOverlaps(gr, hmm)
hits <- data.table(ix = from(hits), state = hmm$state[to(hits)])
hits <- hits[, .N, .(ix, state)]
DT <- dcast(hits, ix ~ state, value.var = "N", fill = 0L)
DT[, ix := NULL]
mcols(gr) <- DT

DT <- as.data.table(gr)
DT[, c("width", "strand") := NULL]
fwrite(DT, file.path(outdir, "wgEncodeBroadHmm_state_count.csv"))

hits <- findOverlaps(clusters, gr)
DT <- cbind(data.table(cluster = as.character(clusters[from(hits)])),
            as.data.table(mcols(gr)[to(hits),]))
DT <- DT[, lapply(.SD, mean), cluster]
DT <- as.matrix(DT, rownames = 1L)
iM <- rowMaxs(DT)
states <- c(colnames(DT), "Mixed")
states <- states[fifelse(iM > 3, apply(DT, 1L, which.max), ncol(DT) + 1L)]
DT <- data.table(cluster = rownames(DT), state = states, key="cluster")
cluster_annot[DT, state := i.state]
# fwrite(DT, file.path(outdir, "meth-cluster_ChromHMM-states.csv"))


# CpG Islands -------------------------------------------------------------

cgi <- AnnotationHub()[["AH5086"]]
cgi <- cgi[seqnames(cgi) %in% paste0("chr", c(1:22, "X", "Y"))]
cgi <- sort(cgi)
cgi$name <- as.character(cgi)
cgi$Relation_to_Island <- "Island"

cgsea <- gaps(cgi)
cgsea <- cgsea[strand(cgsea) == "*" & seqnames(cgsea) %in% paste0("chr", c(1:22, "X", "Y"))]
cgsea <- sort(cgsea)
cgsea$name <- paste0(seqnames(cgsea), ":", start(cgsea), "-", end(cgsea))
cgsea$Relation_to_Island <- "OpenSea"

DT <- mergeByOverlaps(clusters, cgi) |>
    with(data.table(cluster = as.character(clusters),
                    Island = as.character(cgi)))
DT <- DT[, .SD[1], cluster]

DT2 <- mergeByOverlaps(clusters, cgsea)
DT2 <- with(DT2, data.table(cluster = as.character(clusters), OpenSea = as.character(cgsea)))
DT2 <- DT2[, .SD[1], cluster]
DT <- merge(DT, DT2, by = "cluster", all = TRUE)
DT[, Island_Status := fifelse(is.na(Island), "OpenSea", "Island")]
DT[is.na(Island), Island := OpenSea]
DT[, OpenSea := NULL]
setkey(DT, cluster)
cluster_annot <- merge(cluster_annot, DT, all.x = TRUE)
# fwrite(DT, file.path(outdir, "meth-cluster_Islands.csv"))


# TFBS --------------------------------------------------------------------

tfbs <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encRegTfbsClustered/encRegTfbsClusteredWithCells.hg19.bed.gz"
tfbs <- fread(tfbs, sep = "\t", header = FALSE,
              col.names = c("chrom", "start", "end", "name", "score", "sources"))
tfbs <- makeGRangesFromDataFrame(tfbs[score > 600], keep.extra.columns = TRUE)

hits <- findOverlaps(clusters, tfbs, type = "within")
res <- as.data.table(mcols(tfbs)[to(hits), ])
res <- res[, `:=`(cluster = as.character(clusters)[from(hits)],
                  Ncells = stringr::str_count(sources, ",") + 1L,
                  sources = NULL)]
setcolorder(res, "cluster")
fwrite(res, file.path(outdir, "meth-cluster_TFBS.csv"))


# Loops -------------------------------------------------------------------

DT <- fread("data/experiments/Bhattacharyya_loops.csv.gz")#[Source == "H3K27ac"]
loop_mat <- split(DT, by = c("Source", "Cell")) |>
    lapply(function(D) GInteractions(D[, GRanges(seqnames = Chromosome, ranges = IRanges(start = Start1, end = End1))],
                                     D[, GRanges(seqnames = Chromosome, ranges = IRanges(start = Start2, end = End2))])) |>
    vapply(countOverlaps, rep(0L, length(clusters)), query = clusters) |>
    as.data.table()
loop_mat[, cluster := as.character(clusters)]
setkey(loop_mat, cluster)
setcolorder(loop_mat)
fwrite(loop_mat, file.path(outdir, "meth-cluster_loops.csv"))

# Open Chromatin ----------------------------------------------------------

k562_dnase <- "https://www.encodeproject.org/files/ENCFF621ZJY/@@download/ENCFF621ZJY.bed.gz"
k562_dnase <- fread(k562_dnase, select = 1:3, col.names = c("seqnames", "start", "end"))
k562_dnase <- makeGRangesFromDataFrame(k562_dnase, starts.in.df.are.0based = TRUE)
DT <- data.table(cluster = as.character(clusters), K562_DNase = countOverlaps(clusters, k562_dnase), key = "cluster")
cluster_annot[DT, K562_DNase := i.K562_DNase]
# fwrite(DT, file.path(outdir, "meth-cluster_DNase.csv"))

fwrite(cluster_annot, file.path(outdir, "meth-cluster_annotated.csv"))
