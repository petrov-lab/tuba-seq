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
    
### preprocess.py

Filters raw FASTQ files from deep sequencing using a template of the intended barcode setup. Outputs training files for *DADA2* estimation of error rates and clutering files for *DADA2* barcode clustering. Uses nwalign.pyx/nwalign.c for [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) based alignment of reads to the barcode template. 

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