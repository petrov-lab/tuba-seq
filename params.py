

############################## Amplicon Format ################################
# 
# This pipeline expects DNA sequences containing a non-degenerate head, 
# followed by a sequence identifying the sgRNA (enumerated in sgRNA_info.csv),  
# followed by a random nucleotide barcode, followed by a non-degenerate tail.  
#
###############################################################################

head = 'NCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGA'      #  49 
sgID = '........'                                               #   8 
random_barcode = 'AA.....TT.....AA.....'                        #  21
tail = 'ATGCCCAAGAAGAAGAGGAAGG'                                 #  22
                                                                #  --
master_read = head + sgID + random_barcode + tail               # 100


alignment_flank = 22
training_flank = 16             # Number of flanking nucleotides to use for error estimation
cluster_flank = 6               # Number of flanking nucleotides for clustering
min_align_score_frac = 0.5      # lost N cassettes are at ~0.6 - 0.67
allowable_deviation = 4         # Flexibility in length of degenerate region
maxEE = 2                       # Maximum # of Expected Errors

preprocessed_dir = 'preprocessed/'
training_dir = 'training/'
original_dir = 'original/'
fastq_handle = '.fastq'

# Summary


import numpy as np

# Dictionary of reduction operators to use for aggregating DADA2 clusters that are deemed identical (n10 = n0 + n1): 
merge_rules = dict(abundance=np.sum, n0=np.sum, n1=np.sum, nunq=np.sum, pval=np.prod, birth_pval=np.prod, birth_ham=np.min)

sgID_length = len(sgID)
random_barcode_length = len(random_barcode)
barcode_length = sgID_length + random_barcode_length

sample_meta_data_file = ""

#### Panda-seq Parameters   ############
panda_seq = dict(
    merge_dir = 'merged/',
    threshold = 0.6,               # Alignments lower than this value will be discarded (default: 0.6)
    forward_primer = head[1:26],
    reverse_primer = 'AGGTTCTTGCGAACCTCATCACTCGTTGCATCGACCGGTAATGCAGGCAAATTTTGGTGTACGGTCAGTAAATTGGACACCTTCCTCTTCTTCTTGGGCAT'[:56],
    min_overlap = barcode_length,
    min_length = barcode_length + 2*alignment_flank)
