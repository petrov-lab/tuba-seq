#!/usr/bin/env Rscript
#
# DADA2_clustering.R derep/INPUT_FILE.RData 
#
# Clusters a single derep'ed input sample using DADA2. If no input file is 
# provided, this script will assume it is being executed as part of a job array
# and will search the environment variables (verbosely) for a task ID that will
# assign it a job within the `in.directory` folder. This script can take
# several days, sometimes longer, to complete per sample and it is recommended 
# to be used as part of a job array on a high-performance computing cluster.
# Outputs a CSV-format data table entitled INPUT_FILE.csv.gz into the
# `out.directory` folder.
#
#
# A few notes on performance:
#
# 1. Decreasing omega (i.e. a more negative exponent) increases speed. However
#    doing so may over-cluster barcodes. Conversely, a too large omega can be 
#    fixed ex post facto, so don't decrease this value unnecessarily. 
#
# 2. The other parameters have been somewhat optimized for speed on Xeon 
#    architectures, I wouldn't bother tinkering, unless runtimes often 
#    exceed a week.
#
# 3. Multi-threading decreases runtime to ~20-40% in my experience. I recommend
#    broadcasting sample clustering onto a distributed computing system where
#    each sample is clustered by an entire, dedicated compute node.

############################### I/O Options ###################################
in.directory <- "preprocessed"
out.directory <- "clustered"
training.filename <- "trainingDADA.rds"
overwrite.existing.file <- FALSE

############## IMPORTANT CLUSTERING DADA2 OPTIONS #############################
library(dada2, quietly=TRUE)
setDadaOpt( OMEGA_A     = 1e-10,
            USE_QUALS   = TRUE )         # Will ignore Phred Scores if FALSE

########################## PERFORMANCE DADA2 OPTIONS ##########################
setDadaOpt( BAND_SIZE=4,        # Lower values improve performance. This value 
                                # can be as small as `allowable_deviation` 
                                # without compromising alignment quality.
    VECTORIZED_ALIGNMENT=FALSE, # Unproductive for small BAND_SIZE. 
	USE_KMERS=TRUE)             # Default           

####################### I/O SETUP #############################################

# Get file to cluster; use trailing arguments, then try job array IDs
args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 1) {
    derep.file <- args[1];
} else {
    files <- list.files(in.directory);
    dereped.files <- sort(files[grepl(".rds", files)]);
    if (length(dereped.files) == 0)
        stop("No input file provided and no input files in ", in.directory, call. = TRUE)
# Various potential grid-engine task ids:
    slurm.id <- Sys.getenv("SLURM_ARRAY_TASK_ID");
    sge.id <- Sys.getenv("SGE_TASK_ID");
    if (slurm.id != "") {
        derep.base <- dereped.files[[strtoi(slurm.id)]];
        message("Running slurm job #", slurm.id, " of ", length(dereped.files), ".")
    } else if (sge.id != "") {
        derep.base <- dereped.files[[strtoi(sge.id)]];
        message("Running SGE job #", sge.id, " of ", length(dereped.files), ".")
    } else {
        stop("Found input files, but no slurm nor SGE task ID.")
    }
    derep.file <- file.path(in.directory, derep.base);
}
derep <- readRDS(derep.file)

# Load Error Model
if (grepl(".rds", training.filename)) {
  error.model.dadas <- readRDS(training.filename)
  error <- error.model.dadas[[1]]$err_out 
} else if (grepl('.csv', training.filename)) {
  trans <- read.csv(file=training.filename, header=TRUE)
  colnames(trans) <- as.integer(substr(names(trans), start=2, stop=99))
  error <- loessErrfun(as.matrix(trans))
} else {
  stop("training.filename must be a csv or rds file.")
}


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
message("Successfully clustered ", sample.name, " in ", round(proc.time()[['elapsed']]), " seconds.")

###############################################################################
