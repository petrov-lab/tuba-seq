#!/bin/sh

python3 ../preprocess_unadulterated.py

gunzip training/*.gz
gunzip preprocessed/*.gz

../DADA2_derep.R ./
../DADA2_error_training.R training 2> error_training_output_unadulterated.txt
