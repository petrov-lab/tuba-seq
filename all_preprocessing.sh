#!/bin/sh
set -e 

TUBA_SEQ_DIR="$HOME/tuba_seq/"
DATA_DIR=`pwd` 

echo "Analyzing files in $DATA_DIR."

$TUBA_SEQ_DIR/preprocess_unadulterated.py -p $DATA_DIR;
$TUBA_SEQ_DIR/DADA2_derep.R $DATA_DIR;
$TUBA_SEQ_DIR/DADA2_error_training.R $DATA_DIR/training_derep.RData;

echo "All preprocessing completed. If the error training output stats and the fraction of cluster-able reads look good, then you are ready for clustering.";
echo "This task ought to be broadcasted onto a cluster.";
echo "See DADA2_clustering.R for details."

