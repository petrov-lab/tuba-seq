import numpy as np
############################## Amplicon Format ################################
# 
# This pipeline expects DNA sequences containing a non-degenerate head, 
# followed by a sequence identifying the sgRNA (enumerated in sgRNA_info.csv),  
# followed by a random nucleotide barcode, followed by a non-degenerate tail.  
#
###############################################################################
                                                            #  Lengths
head = 'GCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGA'      #  49 
sgID = '........'                                               #   8 
random_barcode = 'AA.....TT.....AA.....'                        #  21
tail = 'ATGCCCAAGAAGAAGAGGAAGGTGTCCAATTTACTGACCGTACACCAAAATTTGCCTGCATTACCGGTCGATGCAACGAGTGATGAGGTTCGCAAGAACCT'[:22]
                                                                #  22
                                                                #  --
master_read = head + sgID + random_barcode + tail               # 100

######################## Preprocessing Parameters #############################

min_align_score_frac = 0.6      # lost N cassettes are at ~0.6 - 0.67
allowable_deviation = 4         # Flexibility in length of degenerate region
alignment_flank = 22
training_flank = 18             # Number of flanking nucleotides to use for error estimation
cluster_flank = allowable_deviation # Number of flanking nucleotides for clustering

preprocessed_dir = 'preprocessed/'
training_dir = 'training/'
original_dir = 'original/'
fastq_handle = '.fastq'

# Dictionary of reduction operators to use for aggregating DADA2 clusters that are deemed identical (n10 = n0 + n1): 
merge_rules = dict(abundance=np.sum, n0=np.sum, n1=np.sum, nunq=np.sum, pval=np.prod, birth_pval=np.prod, birth_ham=np.min)

sgID_length = len(sgID)
random_barcode_length = len(random_barcode)
barcode_length = sgID_length + random_barcode_length

sample_meta_data_file = ""

