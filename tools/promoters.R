# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: promoters.R
# time: 2021-04-20
suppressPackageStartupMessages(library(GenomicFeatures))

# parsing args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste(
        "Usage: Rscript", "promoters.R",
        "txdb.sqlite", "promoters.bed"
    ))
}

txdb <- loadDb(args[1])
pr <- suppressWarnings(promoters(txdb, upstream = 2000, downstream = 400))
pr <- trim(pr)
starts <- start(pr) - 1
starts[starts < 0] <- 0
pr <- data.frame(
    seqnames = seqnames(pr), starts = starts, ends = end(pr)
)

write.table(
    pr, args[2],
    quote = FALSE, sep = "\t", na = "", row.names = FALSE,
    col.names = c("#chrom", "start", "end")
)