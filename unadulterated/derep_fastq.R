#!/usr/bin/Rscript

file <- commandArgs(trailingOnly=TRUE)[1]
library(dada2)

derep <- derepFastq(file, verbose=TRUE)
noExt <- strsplit(fastq, ".fastq")[[1]]
saveRDS(preprocess.derep, file=paste0(noExt, ".rds"))

