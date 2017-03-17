
searches = 8
E_value_threshold = 1e-10
short_title = True
multithread = True

Local_BLASTX = True

# BLAST database
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

temp_fasta = 'temp.fasta'
blast_xml = 'blast_results.xml'


import pandas as pd

def BLASTn_search(seq):
    results_handle = NCBIWWW.qblast("blastn", 'nt', seq, descriptions=2) #, hitlist_size=1)
    blast_records = NCBIXML.parse(results_handle)
    blast_record = next(blast_records)              # Preformed only one BLAST search
    description = blast_record.descriptions[0]   # Enforced as only 1 (using hitlist_size=1)
    return description

def sleuth_DNAs(DNAs):
    counts = DNAs.value_counts()
    top_hits = counts.nlargest(searches)
    top_hits.name = 'reads'    
    with open(temp_fasta, 'w') as f:
        f.write('\n'.join(top_hits.index.values) + '\n')

    blastx_cline = NcbiblastnCommandline(query=temp_fasta, db="nt", evalue=E_value_threshold, outfmt=5, out=blast_xml)
    stdout, stderr = blastx_cline()
    with open(blast_xml) as f:
        blast_records = list(NCBIXML.parse(f))

    descriptions_df = pd.DataFrame([record.descriptions[0].__dict__ for record in blast_records], index=top_hits.index)
    output = pd.merge([top_hits.to_frame(), descriptions_df], axis=1)
    if short_title:
        output['title'] = output['title'].apply(lambda s: s.split('|')[-1])
    return output



def sleuth_DNAs(DNAs):
    counts = DNAs.value_counts()
    top_hits = counts.nlargest(searches)
    top_hits.name = 'reads'    
    if multithread:
        from multiprocessing import Pool
        with Pool(processes=min(searches, 999)) as p:
            descriptions = p.map(BLASTn_search, top_hits.index.values)
    else:
        descriptions = map(BLASTn_search, top_hits.index.values)
    description_df = pd.DataFrame([des.__dict__ for des in descriptions], index=top_hits.index)

    merged = pd.merge([top_hits.to_frame(), description_df], axis=1)
    output = merged.query("expect <= {:}".format(E_value_threshold))
    if short_title:
        output['title'] = output['title'].apply(lambda s: s.split('|')[-1])
    return output

_start, _stop, max_score = NW_fit(b_ref, b_ref, training_flank)
_start, _stop, min_score = NW_fit(to_bytes( ref.replace('A', 'T').replace('G', 'C') ), to_bytes( ref.replace('T', 'A').replace('C', 'G') ), training_flank)

