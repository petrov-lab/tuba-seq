# TuBa-seq Analysis Pipeline #

Command-line data analysis tools for tumor size profiling by TUmor-Barcoded deep sequencing. Please refer to our publication in [Nature Methods](http://www.nature.com/nmeth/index.html) for details. 

##OVERVIEW

This pipeline uses the [DADA2](https://github.com/benjjneb/dada2) de-noising and sample inference algorithm to identify unique DNA barcodes. Sequencing error rates are estimated from the non-degenerate regions of DNA barcodes using a novel approach. While *DADA2* is written in R, the pre-processing and post-processing scripts, and utilities are written in Python. 

##INSTALLATION & USAGE

##CONTACT
---------
Feedback is most appreciated. Please contact me at cmcfarl2 [at] stanford [dot] edu with any questions or suggestions. 

WEBSITE
-------
Documentation and source code can be found at 
[https://bitbucket.org/cdmcfarland/tuba-seq]

## FILE MANIFEST
----------------
    
Often two versions of files exist: a base version and a version with an "_unadulterated" suffix. Base versions are recommended for general use. Unadulterated versions fully replicate the findings presented in [Rodgers et al. (2017)](http://www.nature.com/nmeth/index.html). These versions do not differ meaningful, however we want to provide users a way to fully reproduce our published scientific findings, while also developing a generic pipeline that assimilates best practices as they evolve. Unadulterated versions are un-maintained. 

### preprocess.py

Filters raw FASTQ files from deep sequencing using a template of the intended barcode setup. Outputs training files for *DADA2* estimation of error rates and clutering files for *DADA2* barcode clustering. Uses nwalign.pyx/nwalign.c for [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) alignment of reads to the barcode template for filtering & trimming. 


### postprocess.py

Consolidates *DADA2* clustering outputs. Annotates sgRNAs and DNA benchmarks, and isolates the barcode region of reads. Clusters in separate sequencing runs, belonging to the same mouse, and clusters known to derive from the same DNA barcode (by virtue of the DNA template) are consolidated. Summary statistics are provided. 

### final_processing.py

Creates the final estimates of tumor sizes used in this study. These steps are application specific, and we have converged upon best practices for generic tuba-seq experiments. Mice are annotated by their germ-line genotype, their time of sacrifice, and by the viral pool used to initiate tumors. Viral infections not belonging to the intended pool used are removed (these are very rare events). Estimates any bias in PCR amplification that is correlated with the GC content of barcodes and subtracts this bias. 

### What is this repository for? ###

* Quick summary
* Version
* [This Method ](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact