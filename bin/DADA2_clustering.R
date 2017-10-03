#!/usr/bin/env Rscript
#
# DADA2_clustering.R derep/INPUT_FILE.RData 
#
# Clusters a single derep'ed input sample using DADA2. This script can take
# several days, sometimes longer, to complete per sample and it is recommended 
# to be broadcasted onto a high-performance computing cluster (one job per 
# sample.) Outputs a CSV-format data table entitled INPUT_FILE.csv.gz into the
# `out_directory` folder.
#
#
# A few notes on performance:
#
# 1. Decreasing omega (i.e. a more negative exponent) increases speed. However
#    doing so may over-cluster barcodes. Conversely, a too large omega can be 
#    fixed ex post facto, so don't decrease this value unnecessarily. 
#
# 2. The other parameters have been somewhat optimized for speed on common 
#    architectures, I wouldn't bother with tinkering, unless runtimes often 
#    exceed a week.
#
# 3. Multithreading decreases runtime by ~25-50% in my experience. I recommend
#    broadcasting sample clustering onto a distributed computing system where
#    each sample is clustered by a dedicated compute node (with no other
#    concurrent tasks.)

############################### I/O Options ###################################
out.directory <- "clustered"
training.filename <- "trainingDADA.rds"
overwrite.existing.file <- FALSE

############## IMPORTANT CLUSTERING DADA2 OPTIONS #############################
paired.end.reads <- TRUE
library(dada2, quietly=TRUE)

if (paired.end.reads) {
    setDadaOpt( OMEGA_A     = 1e-2,
                USE_QUALS   = TRUE)         # Will ignore Phred Scores if TRUE
} else {
    setDadaOpt( OMEGA_A     = 1e-10,
                MIN_HAMMING = 2,            # Possibly unnecessary
                MIN_FOLD    = 2,            # Possibly unnecessary
                USE_QUALS   = TRUE)         # Will ignore Phred Scores if TRUE
}
########################## PERFORMANCE DADA2 OPTIONS ##########################
setDadaOpt( BAND_SIZE=4,    # Lower values improve performance. This value can
                            # be as small as `allowable_deviation` in params.py
                            # without compromising alignment quality.
    VECTORIZED_ALIGNMENT=FALSE, # Unproductive for small BAND_SIZE. 
	USE_KMERS=TRUE)             # Default           

####################### I/O SETUP #############################################

args <- commandArgs(trailingOnly=TRUE)
derep.file <- args[1]
derep <- readRDS(derep.file)
error.model.dadas <- readRDS(training.filename)
error <- error.model.dadas[[1]]$err_out 

sample.name <- basename(strsplit(derep.file, ".rds")[[1]][[1]])

dir.create(out.directory, showWarnings=FALSE)
out.file <- file.path(out.directory, paste0(sample.name, ".csv.gz"))
if (file.exists(out.file) && !overwrite.existing.file)
  stop(out.file, " already created, exiting now.", call. = TRUE)

########################## DADA2 Clustering ###################################

message("Clustering ", sample.name, ' to ', out.file, ' ...')
output <- dada( derep, 
                errorEstimationFunction=loessErrfun,  # Avoids useless warning
		        err=error,
		        selfConsist=FALSE,
                multithread=TRUE)       # See Note 3.    

File <- gzfile(out.file, "w")
write.csv(output$clustering, File, quote=FALSE)
close(File)
message("Successfully clustered ", sample.name, ".")

###############################################################################
