#!/usr/bin/env python3
import sys
import logging
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_pi(fasta_input, tsv_output, log_file=None):
    """
    Calculates the isoelectric point (pI) and GRAVY score for protein sequences.
    
    Args:
        fasta_input (str): Path to the input FASTA file containing protein sequences.
        tsv_output (str): Path to the output TSV file.
        log_file (str, optional): Path to a log file.
    """
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(fh)

    logger.info(f"Starting pI calculation for {fasta_input}")
    results = []
    
    try:
        for record in SeqIO.parse(fasta_input, "fasta"):
            # Clean sequence: remove stop codons (*) and non-standard AA (X, B, Z, J)
            seq_str = str(record.seq).replace("*", "")
            # Keep only standard amino acids
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
                logger.warning(f"Skipping sequence {record.id}: {e}")
                continue
        
        logger.info(f"Processed {len(results)} sequences.")
        df = pd.DataFrame(results)
        df.to_csv(tsv_output, sep="\t", index=False)
        logger.info(f"Results saved to {tsv_output}")

    except Exception as e:
        logger.error(f"Fatal error processing {fasta_input}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if 'snakemake' in locals():
        log_path = snakemake.log[0] if snakemake.log else None
        calculate_pi(snakemake.input[0], snakemake.output[0], log_path)
    elif len(sys.argv) >= 3:
        calculate_pi(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python calculate_PI.py <input.fasta> <output.tsv>")
        sys.exit(1)