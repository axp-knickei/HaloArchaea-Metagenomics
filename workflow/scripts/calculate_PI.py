import sys
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

def calculate_pi(fasta_input, tsv_output):
    results = []
    
    for record in SeqIO.parse(fasta_input, "fasta"):
        # Clean sequence: remove stop codons (*) and non-standard AA (X, B, Z, J)
        seq_str = str(record.seq).replace("*", "")
        clean_seq = "".join([x for x in seq_str if x not in "XBZJ"])
        
        # Skip very short peptides
        if len(clean_seq) < 10:
            continue
            
        try:
            analyser = ProteinAnalysis(clean_seq)
            pI = analyser.isoelectric_point()
            gravy = analyser.gravy()
            mw = analyser.molecular_weight()
            
            results.append({
                "protein_id": record.id,
                "pI": pI,
                "gravy": gravy,
                "molecular_weight": mw,
                "length": len(clean_seq)
            })
        except Exception as e:
            continue
            
    df = pd.DataFrame(results)
    df.to_csv(tsv_output, sep="\t", index=False)

if __name__ == "__main__":
    calculate_pi(sys.argv[1], sys.argv[2])