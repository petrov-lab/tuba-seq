#!/usr/bin/Rscript

file <- commandArgs(trailingOnly=TRUE)[1]
library(dada2)

preprocess.derep <- derepFastq(file, verbose=TRUE)
noExt <- strsplit(file, ".fastq")[[1]][[1]]
message(noExt)
message(paste0(noExt, ".rds"))
saveRDS(preprocess.derep, file=paste0(noExt, ".rds"))

