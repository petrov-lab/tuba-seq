#!/usr/bin/env Rscript
# 
# DADA2_error_training.R derep_training_directory
#
#' Develops a model of read errors for barcode clustering. The training directory
#' contains all the de-repped training files (*.rds) to be fit. Training files 
#' only contain the nucleotides flanking the DNA barcodes, which should all be
#' identical save for read errors. Thus, very conservative, iterative clustering 
#' is used to ascertain error rates. Outputs the error model, a graph of its 
#' rates, and summary measures of the error rates of the Illumina run. 

library(dada2, quietly=TRUE)
library(tools)

args <- commandArgs(trailingOnly=FALSE)
script.filename <- gsub('--file=', "", grep('--file=', args, value=TRUE)[[1]])
script.directory <- dirname(script.filename)
derep.directory <- commandArgs(trailingOnly=TRUE)[1] 

omega <- 1e-300         # 'The alpha & the omega' parameter of barcode clustering
# Omega is very small because multiple training clusters are highly unlikely.
# This value (1e-300) seems to minimize bias better than the optimal clustering 
# value (omega=1e-10) and no clustering (omega=0) presumably because there are 
# very-rare mutations/contaminations that should fall into separate clusters. 
paired.end.reads <- TRUE

######################## I/O Options ########################################## 

training.filename <- "trainingDADA.rds"
error.graph.filename <- "ErrorModel.pdf"

############## IMPORTANT TRAINING DADA2 OPTIONS ##############################
setDadaOpt( 
            MIN_HAMMING=2,  # Possibly unnecessary; 2 was used in Rodgers et al
            MIN_FOLD=2,
            USE_QUALS=TRUE) # Will ignore Phred Scores if TRUE

########################## PERFORMANCE DADA2 OPTIONS ##########################
setDadaOpt(
  BAND_SIZE=4,  # Lower values improve performance. This value can be as small
                # as `allowable_deviation` without compromising alignments.
  VECTORIZED_ALIGNMENT=FALSE, # Unproductive for small BAND_SIZE. 
  USE_KMERS=TRUE,       # Default           
  MAX_CONSIST=200)      # Convergence can take many iterations, doesn't matter.       

####################### Load Files ############################################

files <- list.files(derep.directory)
dereped.files <- sort(files[grepl(".rds", files)])
dereps <- lapply(file.path(derep.directory, dereped.files), readRDS)
sample.filenames <- sapply(dereped.files, basename)
names(dereps) <- sapply(sample.filenames, file_path_sans_ext)

####################### ERROR TRAINING ########################################
message("Running Error-training DADA2 on ", length(dereps), " samples with omega = ", omega)

dadas <- dada(  dereps,
                OMEGA_A=omega,
                selfConsist=TRUE,
                err=NULL,
                MAX_CONSIST=200,
                multithread=TRUE)

message("Saving model to ", training.filename)
saveRDS(dadas, file=training.filename)

####################### OUTPUT SUMMARY ########################################

dfs <- lapply(dadas, function(cluster) cluster$clustering)

message("Graphing the LOWESS-regressed error model. Generally observed Phred error rates exceed expectation (red) by 2-3 fold.")
pdf(error.graph.filename)
plotErrors(dadas[[1]], nominalQ=TRUE)

largest.cluster.fractions <- sapply(dfs, function(df) df$abundance[[1]]/sum(df$abundance) )

message("")
message("Error training Summary for omega = ", omega)
message("The mean size of the largest cluster is ", round(mean(largest.cluster.fractions)*1e2, 1), "%.")
message("This cluster should dominate (i.e. exceed 95% of the total) as, ideally, there is only 1 training cluster.")

trans <- lapply(dadas, function(df) df$trans)
all_trans <- Reduce("+", trans)
model <- dadas[[1]]$err_out 
correct_nuc <- c('A2A', 'C2C', 'G2G', 'T2T')
corrects <- colSums(all_trans[correct_nuc,])
totals <- colSums(all_trans)

model.errs <- colMeans(model[correct_nuc,])

message("Estimated a total error rate of ", sprintf("%.4f%%",(sum(totals) - sum(corrects))*100/sum(totals) ), " (0.1 - 0.5% is typical.) By Phred Score:")

#cohorts <- if(paired.end.reads) c(0, 20, 40, 50, 56, 58, 60, 61, 62, 63) else c(0, 20, 30, 36, 38, 39, 40, 41)
cohorts <- c(0, 20, 30, 36, 38, 39, 40, 41)

message("Phred Range | Error Rate | LOWESS-Average | Millions of Nucleotides in Range")
for (i in 1:(length(cohorts)-1)) {
  Total <- sum(totals[(cohorts[[i]]+1):cohorts[[i+1]]]);
  Correct <- sum(corrects[(cohorts[[i]]+1):cohorts[[i+1]]]);
  
  err.rate <- 1 - mean(model.errs[(cohorts[[i]]+1):cohorts[[i+1]]]);
  message(sprintf("%d-%d           %.3f%%       %.3f%%       %.1f", 
        cohorts[[i]], cohorts[[i+1]] - 1, (Total - Correct)*100/Total, err.rate*100, Total*1e-6));
}

