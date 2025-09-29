#!/usr/bin/env python3
"""
find_mutation_from_peptide_to_excel.py

Given a single mutant peptide, scan wild-type FASTA sequences
to find where it aligns except for one (or more) mismatches.
Reports the location and the mutation(s), saves results to Excel.
"""

from Bio import SeqIO
import pandas as pd

# === User parameters ===
fasta_file = "/Users/ayalinchjonathan/.Trash/uniprotkb_proteome_UP000005640_2025_09_29.fasta"
mutant_peptide = "ILDTAGKEEY"   # <-- hardcoded peptide
output_file = "mutation_scan_results.xlsx"
max_mismatches = 3  # allow up to N mismatches

# === Scan all proteins ===
results = []

for record in SeqIO.parse(fasta_file, "fasta"):
    protein_seq = str(record.seq)
    pep_len = len(mutant_peptide)

    for i in range(len(protein_seq) - pep_len + 1):
        window = protein_seq[i:i+pep_len]

        mismatches = [
            (pos, window[pos], mutant_peptide[pos])
            for pos in range(pep_len)
            if window[pos] != mutant_peptide[pos]
        ]

        if 0 < len(mismatches) <= max_mismatches:
            results.append({
                "Peptide": mutant_peptide,
                "Protein_ID": record.id,
                "Protein_Description": record.description,
                "Match_Start": i+1,
                "Match_End": i+pep_len,
                "Wildtype_Window": window,
                "Num_Mismatches": len(mismatches),
                "Mutations": "; ".join(
                    [f"{ref}{i+pos+1}{mut}" for pos, ref, mut in mismatches]
                )
            })

# === Save to Excel ===
if results:
    df = pd.DataFrame(results)
    df.to_excel(output_file, index=False)
    print(f"Saved {len(df)} matches to {output_file}")
else:
    print("No matches found.")