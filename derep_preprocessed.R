#!/usr/bin/Rscript

path <- paste0(commandArgs(trailingOnly=TRUE)[1], '/')
#path <- './'
library(dada2, quietly=TRUE) 
compress = 'xz'

preprocessed_path <- paste0(path, 'preprocessed/')

fns <- list.files(preprocessed_path)
preprocessed_fastqs <- sort(fns[grepl(".fastq", fns)])
    
message("\nNow Preprocessed FASTQs.\n")

out_path <- paste0(path, "derep/")
dir.create(out_path, showWarnings=FALSE)

for (fastq in preprocessed_fastqs) {
  derep <- derepFastq(paste0(preprocessed_path, fastq), verbose=TRUE)
  shortName <- basename(strsplit(fastq, ".fastq")[[1]][[1]])
  save(derep, file=paste0(out_path, shortName, ".RData"), compress=compress)
}

unlink(preprocessed_path, recursive=TRUE)

