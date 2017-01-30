#!/usr/bin/Rscript

path <- paste0(commandArgs(trailingOnly=TRUE)[1], '/')
path <- './'
library(dada2, quietly=TRUE) 

############### Input parameters ##############################################

compression <- 'xz'
training_path <- paste0(path, 'training/')
preprocessed_path <- paste0(path, 'preprocessed/')
remove_fastq_files <- TRUE
error_training_outfile <- paste0(path, "training_derep.RData")
clustering_out_path <- paste0(path, "derep/")

###############################################################################

fns <- list.files(training_path)
training_fastqs <- sort(fns[grepl(".fastq", fns)])

fns <- list.files(preprocessed_path)
preprocessed_fastqs <- sort(fns[grepl(".fastq", fns)])
    
message("Found ", length(training_fastqs), " error training files.")
message("Found ", length(preprocessed_fastqs), " clustering files.") 

sample_file_names <- sapply(strsplit(training_fastqs, "/"), tail, n=1)
sample_names <- sapply(strsplit(sample_file_names, ".fastq"), `[`, 1)

# No filtering 
# preprocessing.py ensures that maxN=0, maxEE=desired & that reads are trimmed
derep <- lapply(paste0(training_path, training_fastqs), derepFastq, verbose=TRUE)
names(derep) <- sample_names

save(derep, file=error_training_outfile, compress=compression)

message("\nNow dereping Clustering FASTQs.\n")

dir.create(clustering_out_path, showWarnings=FALSE)

for (fastq in preprocessed_fastqs) {
  derep <- derepFastq(paste0(preprocessed_path, fastq), verbose=TRUE)
  shortName <- basename(strsplit(fastq, ".fastq")[[1]][[1]])
  save(derep, file=paste0(out_path, shortName, ".RData"), compress=compression)
}

if (remove_fastq_files) {
  message("Removing the preprocessed FASTQ files.")
  unlink(training_path, recursive=TRUE)
  unlink(preprocessed_path, recursive=TRUE)
}
