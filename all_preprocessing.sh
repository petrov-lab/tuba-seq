#!/bin/sh
set -e 

TUBA_SEQ_DIR="~/tuba-seq/"
DATA_DIR=pwd; 


$TUBA_SEQ_DIR/preprocessing.py -p $DATA_DIR;
$TUBA_SEQ_DIR/DADA2_derep.R $DATA_DIR;
$TUBA_SEQ_DIR/DADA2_error_training.R $DATA_DIR;

echo "All preprocessing completed. If the error training output stats and the fraction of cluster-able reads look good, then you are ready for clustering.";
echo "This task ought to be broadcasted onto a cluster.";
echo "Generally, each job in the job array takes a form of: DADA2_clustering.R derep/EACH_SAMPLE_FILE.RData -m error_model.RData"
echo "...but read the DADA2_clustering.R usage "


