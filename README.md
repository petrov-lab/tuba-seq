# TuBa-seq Analysis Pipeline #

Command-line data analysis tools for tumor size profiling by **Tu**mor-**Ba**rcoded deep sequencing. Please refer to our publication in [Nature Methods](http://www.nature.com/nmeth/index.html) for details. 

##OVERVIEW

This pipeline uses the [DADA2](https://github.com/benjjneb/dada2) de-noising and sample inference algorithm to identify unique DNA barcodes. Sequencing error rates are estimated from the non-degenerate regions of DNA barcodes using a novel approach. 

![Tuba-seq Analysis Pipeline.png](https://bitbucket.org/repo/Mjxqa5/images/12810822-Tuba-seq%20Analysis%20Pipeline.png)
##INSTALLATION & USAGE

While *DADA2* is written in R, the pre-processing and post-processing scripts, and utilities are written in Python. Both require Python 3.2 and R 3.2 or later (what a fun coincidence). _DADA2_ should install from **Bioconductor** as follows: 

```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
## Other R packages used for summarizing output to install:
biocLite("shortRead")
biocLite("ggplot2")
```
See [DADA2 Installation](https://benjjneb.github.io/dada2/dada-installation.html) for details. 

The python scripts leverage 

##CONTACT
---------
Feedback is most appreciated. Please contact me at cmcfarl2 [at] stanford [dot] edu with any questions or suggestions. 

WEBSITE
-------
Documentation and source code can be found at https://bitbucket.org/cdmcfarland/tuba-seq

## FILE MANIFEST
----------------
    
Often two versions of files exist: a base version and a version with an "_unadulterated" suffix. Base versions are recommended for general use, while unadulterated versions replicate results published in [Rodgers et al. (2017)](http://www.nature.com/nmeth/index.html). Behavior of these two versions should remain similar, however we want to provide users a way to reproduce our published findings, while also developing a more generic pipeline that incorporates evolving best practices. Unadulterated versions are un-maintained. 

### preprocess.py

Filters raw FASTQ files from deep sequencing using a template of the intended barcode setup. Outputs training files for *DADA2* estimation of error rates and clutering files for *DADA2* barcode clustering. Uses nwalign.pyx/nwalign.c for [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) alignment of reads to the barcode template for filtering & trimming. 


### DADA2_derep.R

De-replicates training & clustering files for *DADA2*. Separating this first step from the algorithm saves space and time (when runs are repeated). 

### DADA2_error_training.R

Estimates a model of sequencing errors from the 

### postprocess.py

Consolidates *DADA2* clustering outputs. Annotates sgRNAs and DNA benchmarks, and isolates the barcode region of reads. Clusters in separate sequencing runs, belonging to the same mouse, and clusters known to derive from the same DNA barcode (by virtue of the DNA template) are consolidated. Summary statistics are provided. 

### final_processing.py

Creates the final estimates of tumor sizes used in this study. These steps are application specific, and we have converged upon best practices for generic tuba-seq experiments. Mice are annotated by their germ-line genotype, their time of sacrifice, and by the viral pool used to initiate tumors. Viral infections not belonging to the intended pool used are removed (these are very rare events). Estimates any bias in PCR amplification that is correlated with the GC content of barcodes and subtracts this bias.