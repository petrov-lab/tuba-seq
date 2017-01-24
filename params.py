
def get_master_read():
    head = 'NCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGA'  #  49
    sgID = '........'                                           #   8 
    barcode = 'AA.....TT.....AA.....'                           #  21
    tail = 'ATGCCCAAGAAGAAGAGGAAGG'                             #  22
                                                                # ---
    master_read = head + sgID + barcode + tail                  # 100
    return master_read

master_read = get_master_read()
training_flank = 16             # Number of flanking nucleotides to use for error estimation
cluster_flank = 6               # Number of flanking nucleotides for clustering
min_align_score_frac = 0.6      # lost N cassettes are at ~0.6 - 0.67
allowable_deviation = 4         # Flexibility in length of degenerate region
maxEE = 2                       # Maximum # of Expected Errors

preprocessed = 'preprocessed/'
training = 'training/'
fastq_handle = '.fastq'

# Summary

hist_number = 3

