#!/usr/bin/env Rscript
doc <- "
Compute CGI gene overlaps and distance

Usage: gene_cgi_overlap.R [--radius=<bp>] <genes> <output>

Options:
  --radius=<bp>   Distance around tss to check [default: 1000000]
"
library(docopt)
argv <- docopt(doc)

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicAlignments)
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(AnnotationHub)
})

hg19 <- BSgenome.Hsapiens.UCSC.hg19
chr.sizes <- seqlengths(hg19)[1:24]
radius <- as.numeric(argv$radius)

genes <- fread(argv$genes)
genes[, Relation_to_Island := "Gene"]
setnames(genes, "id", "name")
genes <- sort(makeGRangesFromDataFrame(genes[chrom %in% names(chr.sizes)], keep.extra.columns = TRUE))
genome(genes) <- "hg19"
seqlengths(genes) <- chr.sizes[seqlevels(genes)]
tss <- resize(genes, 1L, fix = "start")

cgi <- AnnotationHub()[["AH5086"]]
cgi <- sort(cgi[seqnames(cgi) %in% paste0("chr", c(1:22, "X", "Y"))])
cgi$name <- as.character(cgi)

cgsea <- gaps(cgi)
cgsea <- cgsea[strand(cgsea) == "*"]
cgsea$name <- paste0(seqnames(cgsea), ":", start(cgsea), "-", end(cgsea))
cgsea$Relation_to_Island <- "OpenSea"
cgi$Relation_to_Island <- "Island"
world <- sort(c(cgi, cgsea))

hits <- findOverlaps(world, genes)
overlap_dt <- data.table(gene = genes$name[to(hits)], island = world$name[from(hits)])

hits <- findOverlaps(world, trim(tss + 1e6), type = "within")
hits_dt <- data.table(gene = genes$name[to(hits)],
                      island = world$name[from(hits)],
                      Relation_to_Island = world$Relation_to_Island[from(hits)],
                      gene_width = width(genes)[to(hits)],
                      gene_dir = fifelse(as.vector(strand(genes))[to(hits)] == "+", 1L, -1L),
                      dist_start = start(world)[from(hits)] - start(tss)[to(hits)],
                      dist_end = end(world)[from(hits)] - start(tss)[to(hits)],
                      overlap = FALSE)
hits_dt[overlap_dt, on = c("gene", "island"), overlap := TRUE]
hits_dt[, `:=`(dist_start = gene_dir * dist_start, dist_end = dist_end * gene_dir) ][
        , upstream := pmin(dist_start, dist_end) < 0 ][
        , dist := fcase(!upstream & overlap, 0L,
                         upstream & overlap, pmin(dist_start, dist_end),
                         upstream          , pmax(dist_start, dist_end),
                        !upstream          , pmin(dist_start, dist_end) - gene_width)][
        , c("gene_dir", "gene_width") := NULL ]
fwrite(hits_dt, argv$output)
