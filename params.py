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

###################### Paired-end Merging Parameters ##########################
primer_length = 18

ASCII_BASE = 33
def int_to_phred(x): 
    return chr(x+ASCII_BASE)

max_PHRED = 94
merge_params = dict(
    merge_dir = 'merged/',
    threshold = 0.6,                # Alignments lower than this value will be discarded 
    min_overlap = barcode_length,
    min_length = barcode_length + 2*alignment_flank,
    max_phred = int_to_phred(max_PHRED)    # MAX_PHRED for DADA2 is 62, i.e. '_' 
)

# When merging paired-end reads, a combined PHRED score must be assigned to each
# nucleotide. The various merging algorithms choose different methods to combine
# scores. From the trained error model, we can adjudicate whether the combined 
# PHRED score faithfully describes the true error rate. For example, a combined
# PHRED score of 40 *ought* to yield the expected base 99.99% of the time. 
# Panda-seq appears to be too conservative, i.e. a PHRED score of 40 is correct 
# >99.99% of the time. Other merging algorithms, including the recommended 
# algorithm PEAR, simply add PHRED scores and appear to be too aggressive, i.e.
# a PHRED score of 40 is correct <99.99% of the time. However, overly-aggressive
# approaches can be rectified ex post facto. This can be done, crudely, with the
# mapping below: PHRED scores 1-60 are divided by three, while PHRED scores 
# above 60 are not compressed. 

# DADA2 is designed to automatically determine the true error rate of each PHRED
# score, however it makes a linear assumption during read-dereping that is 
# inappropriate when the true error rate and reported error rates disagree by 
# several decades--as they often will without the re-mapping below.  

compress_PHRED = True
attenuation_rate = 3
attenuation_cap = 60
max_attenuated = int(attenuation_cap/attenuation_rate)
PHRED_compressor = {ASCII_BASE + i:(ASCII_BASE + int(i/attenuation_rate)) if i < attenuation_cap else (ASCII_BASE + i - attenuation_cap + max_attenuated) for i in range(max_PHRED)}


