# TuBa-seq Analysis Pipeline #

Command-line tumor-calling & data analysis tools for tumor size profiling by **Tu**mor-**Ba**rcoded ultra-deep sequencing. Please refer to our publication in [Nature Methods][1] for details. 

## OVERVIEW

This package includes a pipeline for processing Illumina &reg Single-end or Paired-end reads (Paired-end is recommended) into unique barcode pileups and their size in Absolute Cell Number, as described in our [original Tuba-Seq paper][1]. Additional tools are provided, including: 

1. Identification of contaminating DNA and a BLAST-search based investigation of the source of contaminating DNA;
2. Analysis of the sequencing error-rate of the Illumina run;
3. Analysis of barcode randomness, total diversity of the viral library, and evidence for cross-contamination within the library;
4. Statistical tools to summarize tumor size distributions and mice cohorts--including identification of outliers;
5. Graphing functions to illustrate size distributions.

This pipeline uses the [DADA2][2] de-noising and sample inference algorithm to identify unique DNA barcodes. Because tumors differ in size by >1,000-fold, it is difficult to delineate small tumors from recurrent read errors ([this issue has been well-described previously][3]). As such, DADA2 statistical modeling is computationally-intesnive and tools for broadcasting this step onto a distributed computing cluster are provided. Sequencing error rates are estimated from the non-degenerate regions of DNA barcodes using a novel approach. The pipeline is also extensible to other barcode clustering algorithms, including direct support for [Bartender][4], which clusters barcodes in a fraction of the time. However, in our experience DADA2 is necessary for faithful tumor calling. 

Proficiency with python, its scientific computing modules, command-line scripting, distributed clusters, and patience is recommended. 

## INSTALLATION

Installation from source is intended for POSIX-stle Operating Systems (e.g. Linux, darwin). The prequisites are:
* Python 3.2 or later
* R 3.2 or later
* gcc or clang

Executable scripts, by virtue of their shebang, assume Python3 & and R are within the user's `$PATH`. The following commands will generally complete installation: 

```
#Download Source Code
git clone https://github.com/petrov-lab/tuba-seq.git  
cd tuba-seq

#Install the tuba-seq Python package 
./setup.py install

#Install DADA2 from source
./INSTALL_DADA2.R

#Add command-line tools to User's Path (optional)
cp bin/* ~/bin
```

To merge paired-end reads, [PEAR][5] is recommended. It can be installled by following instructions [here][5]. This pipeline is in a very nascent state, Troubleshooting should be dealt with by inspection of the install scripts or by contacting me. 

## CONTACT

Feedback is most appreciated. Please contact me at cmcfarl2 [at] stanford [dot] edu with any questions or suggestions. 

## WEBSITE

Documentation and source code can be found at https://github.com/petrov-lab/tuba-seq

## FILE MANIFEST

### bin/ 

Contains all executables. Python script usage can always be accessed via the `--help` command, e.g.
```
bin/preprocess.py --help
```
R scripts are self-documented and several very-simple shell scripts are provided to illustrate the typical pipeline workflow. They're best understood by inspection. 

### unadulterated

Contains the version of this code that was used to generate the data published in [Rodgers et al. (2017)][1]. Behavior of this version remains similar to the modern version, but it is not recommended for use (other than for reproducing our published findings). This version lacks considerable functionality and does not incorporate what we believe to be the, evolving, best-practices.  

### test/

Test suite. 

### tuba_seq/

Directory of shared python modules.

[1]: https://www.nature.com/nmeth/journal/v14/n7/full/nmeth.4297.html "A quantitative and multiplexed approach to uncover the fitness landscape of tumor suppression in vivo"
[2]: https://github.com/benjjneb/dada2 "DADA2 Public Repository"
[3]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0999-4 "Reproducibility of Illumina platform deep sequencing errors allows accurate determination of DNA barcodes in cells" 
[4]: https://www.biorxiv.org/content/early/2016/08/10/068916 "Bartender: an ultrafast and accurate clustering algorithm to count barcode and amplicon reads"
[5]: https://sco.h-its.org/exelixis/web/software/pear/ "PEAR - Paired-End reAd mergeR"
