#!/usr/bin/env Rscript
# 
# INSTALL_DADA2.R 
#
#' Installs DADA2 from source (github repo). Don't install the binary version--you'll need gcc/clang for the python build anyways!

source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
biocLite("ShortRead", suppressUpdates = FALSE)
biocLite("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2")
