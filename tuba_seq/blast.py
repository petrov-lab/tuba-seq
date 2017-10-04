"""Perform NCBI's BLAST-searching on DNA sequences.

Leverages biopython's Blast sub-package, except always queries the entire 'nt'
database and simplifies output.

By default, web-based searching is used. However, if you have a local 'nt'
database and the NCBI BLAST+ suite, local searches are possible. This database
must be setup as described in 
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc89

"""

import pandas as pd
import os

# BLAST database
from Bio.Blast import NCBIWWW, NCBIXML

E_value_threshold = 1e-10

def web_BLASTn_search(seq, null_result=dict(e=pd.np.inf, title='No Match', bits=0, score=0, num_alignments=0, accession='n.a.')):
    """Web-based BLASTn query of nucleotide sequence `seq` against nt database. 
    Returns dictionary only of best match.
    """
    import warnings
    warnings.warn("Web-based BLAST queries take a highly-variable length of time.")
    try:
        results_handle = NCBIWWW.qblast("blastn", 'nt', seq, hitlist_size=1)
    except ValueError as err:
        print("NCBI server rejected query. Will continue...")
        return null_result
    blast_records = NCBIXML.parse(results_handle)
    blast_record = next(blast_records)              # Preformed only one BLAST search
    des = blast_record.descriptions
    return des[0].__dict__ if des else null_result 

def local_BLASTn_searches(seqs, temp_fasta_file='temp.fasta', BLAST_xml_file='blast_results.xml'):
    """BLAST+ based BLASTn query of iterable nucleotide sequences `seqs` against
    nt database. Your local computer must have BLAST+ tools installed and the nt
    database downloaded. Web-based searching avoids this problem and requires no
    configuring, so I wouldn't use this function unless you are performing a lot 
    of searches. This function is a wrapper around Biopython's BLAST tools. 
    Ensure that the Bio.Blast.Applications.NcbiblastnCommandline class is working
    properly before use. See:

http://biopython.org/DIST/docs/api/Bio.Blast.Applications.NcbiblastnCommandline-class.html

    Returns pandas.DataFrame structure of each nucletoide sequence's best match.
"""
    from Bio.Blast.Applications import NcbiblastnCommandline
    pd.Series(seqs).to_csv(temp_fasta_file, index=False)
    blastx_cline = NcbiblastnCommandline(query=temp_fasta_file, db="nt", evalue=E_value_threshold, outfmt=5, out=BLAST_xml_file, num_threads=8)
    stdout, stderr = blastx_cline()
    os.remove(temp_fasta_file)
    with open(BLAST_xml_file) as f:
        blast_records = list(NCBIXML.parse(f))
    return pd.DataFrame([record.descriptions[0].__dict__ for record in blast_records if record.descriptions], index=seqs)

def sleuth_DNAs(DNAs, local_blast=False):
    """Uses NCBI's BLAST tools to sleuth all known nucleotides for the biological
    origin of unexpected sequences. Returns pandas.DataFrame structure of found 
    matches. 
"""
    output = pd.DataFrame(local_BLASTn_searches(DNAs) if local_blast else list(map(web_BLASTn_search, DNAs)), index=DNAs) 
    output = output.query('e <= {:}'.format(E_value_threshold))
    output['short_title'] = output['title'].apply(lambda s: s.split('|')[-1])
    return output

