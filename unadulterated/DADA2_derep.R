!/usr/bin/env Rscript

directory <- paste0(commandArgs(trailingOnly=TRUE)[1], '/')
library(dada2)
library(base)
############### Input parameters ##############################################

training.path <- paste0(directory, 'training/')
preprocessed.path <- paste0(directory, 'preprocessed/')

###############################################################################

fns <- list.files(training.path)
training.fastqs <- sort(fns[grepl(".fastq", fns)])

fns <- list.files(preprocessed.path)
preprocessed.fastqs <- sort(fns[grepl(".fastq", fns)])
    
message("Found ", length(training.fastqs), " error training files.")
message("Found ", length(preprocessed.fastqs), " clustering files.") 

# No filtering 
# preprocessing.py ensures that maxN=0, maxEE=desired & that reads are trimmed

for (fastq in preprocessed.fastqs) {
  preprocess.derep <- derepFastq(paste0(preprocessed.path, fastq), verbose=TRUE)
  shortName <- basename(strsplit(fastq, ".fastq")[[1]][[1]])
  saveRDS(preprocess.derep, file=paste0(preprocessed.path, shortName, ".rds"))
}

for (fastq in training.fastqs) {
  training.derep <- derepFastq(paste0(training.path, fastq), verbose=TRUE)
  shortName <- basename(strsplit(fastq, ".fastq")[[1]][[1]])
  saveRDS(training.derep, file=paste0(training.path, shortName, ".rds"))
}

