#!/usr/bin/Rscript
# 
# DADA2_error_training.R DEREP_TRAINING_SEQUENCES.RData 
#
#' The original data set contains very long column headers. This function
#' does a keyword search over the headers to find those column headers that
#' match a particular keyword, e.g., mean, median, etc.
#' @param x The data we are querying (data.frame)
#' @param v The keyword we are searching for (character)
#' @param ignorecase Should case be ignored (logical)
#' @return A vector of column names matching the keyword 

args <- commandArgs(trailingOnly=FALSE)
script.filename <- gsub('--file=', "", grep('--file=', args, value=TRUE)[[1]])
script.directory <- dirname(script_filename)
derep.directory <- commandArgs(trailingOnly=TRUE)[1] 
working.directory <- dirname(derep_filename)

######################## I/O Options ########################################## 

training.filename <- "trainingDADA.RData"
error.graph.filename <- "ErrorModel.pdf"

############## IMPORTANT TRAINING DADA2 OPTIONS ##############################
setDadaOpt( OMEGA_A=1e-300,
# Omega is very small because multiple training clusters are highly unlikely.
# This value (1e-300) seems to minimize bias better than the optimal clustering 
# value (omega=1e-10) and no clustering (omega=0) presumably because there are 
# very-rare mutations/contaminations that should fall into separate clusters. 
            MIN_HAMMING=2,  # Possibly unnecessary; 2 was used in Rodgers et al
            MIN_FOLD=2,
            USE_QUALS=TRUE) # Will ignore Phred Scores if TRUE

########################## PERFORMANCE DADA2 OPTIONS ##########################
setDadaOpt( BAND_SIZE=4,    # Lower values improve performance. This value can
                            # be as small as `allowable_deviation` in params.py
                            # without compromising alignment quality.
    VECTORIZED_ALIGNMENT=FALSE, # Unproductive for small BAND_SIZE. 
	USE_KMERS=TRUE,         # Default           
    multithread=TRUE,
    MAX_CONSIST=200)        # Convergence can take a while--doesn't matter.       

####################### Load Files ############################################

library(dada2, quietly=TRUE)
library(tools)
#load(derep_filename)

files <- list.files(derep.directory)
dereped.files <- sort(fns[grepl(".rds", fns)])
dereps <- lapply(file.path(derep.directory, dereped.files), readRDS)
sample.filenames <- sapply(basename, dereped.files)
names(dereps) <- sapply(file_path_sans_ext, sample.filenames)

####################### ERROR TRAINING ########################################

dadas <- dada(  dereps, 
                selfConsist=TRUE,
                err=NULL,
                MAX_CONSIST=200)

trained.error <- dadas[[1]]$err_out
saveRDS(trained.error, file=file.path(working.directory, training.filename))

####################### OUTPUT SUMMARY ########################################
#library(ggplot2, quietly=TRUE)
#source(paste0(script_directory, "R/multiplot.R"))

dfs <- lapply(dadas, function(cluster) cluster$clustering)

message("Graphing the LOWESS-regressed error model. Generally observed Phred error rates exceed expectation (red) by 2-3 fold.")
pdf(file.path(working.directory, error.graph.filename))
plotErrors(dadas[[1]], nominalQ=TRUE)

largest.cluster.fractions <- sapply(dfs, function(df) df$abundance[[1]]/sum(df$abundance) )

message("The mean size of the largest cluster is ", round(mean(largest.cluster.fraction)*1e2, 1), "%.")
message("This cluster should dominate (i.e. exceed 95% of the total) as, ideally, there is only 1 training cluster.")

correct <- sum(sapply(dadas, function(df) as.numeric(rowSums(df[['trans']])[c("A2A", "T2T", "C2C", "G2G")])))
all <- sum(sapply(dadas, function(df) as.numeric(sum(df$trans))))
message("Estimated an error rate of ", round((all - correct)/all*1e2, 3), "% (0.1 - 0.5% is typical.)")

