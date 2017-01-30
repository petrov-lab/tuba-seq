

############################## Amplicon Format ################################
# 
# This pipeline expects DNA sequences containing a non-degenerate head, 
# followed by a sequence identifying the sgRNA (enumerated in sgRNA_info.csv),  
# followed by a random nucleotide barcode, followed by a non-degenerate tail.  
#
###############################################################################

head = 'NCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGA'  #  49
sgID = '........'                                           #   8 
barcode = 'AA.....TT.....AA.....'                           #  21
tail = 'ATGCCCAAGAAGAAGAGGAAGG'                             #  22
                                                                # ---
master_read = head + sgID + barcode + tail                  # 100

training_flank = 16             # Number of flanking nucleotides to use for error estimation
cluster_flank = 6               # Number of flanking nucleotides for clustering
min_align_score_frac = 0.6      # lost N cassettes are at ~0.6 - 0.67
allowable_deviation = 4         # Flexibility in length of degenerate region
maxEE = 2                       # Maximum # of Expected Errors

preprocessed = 'preprocessed/'
training = 'training/'
fastq_handle = '.fastq'

# Summary

hist_number = 30

