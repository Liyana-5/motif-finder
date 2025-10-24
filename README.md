## Motif Finder in FASTA Files

This Python script (find_motif_fasta_txt.py) searches for specific amino acid motifs in protein sequences from a FASTA file and extracts relevant information for downstream analysis.

# Features

- Reads a FASTA file containing protein sequences.
- Reads a tab-delimited product description file with gene IDs and descriptions.
- Searches for a defined amino acid motif in each sequence.
- Extracts the gene ID, product description, and the last 30 amino acids of sequences containing the motif.
- Highlights the motif in lowercase in the last 30 amino acids.
- Outputs results to a tab-delimited file.

# Usage 
**python finding_motif_fasta_txt.py <fasta_file> <description_file> <output_file>**
Example files used in this script is given in the repository as:
 - Fasta file : Tgondii_Proteins.fasta
 - txt file : Tgondii_product_descriptions.txt**

**These files are taken from:https://veupathdb.org/veupathdb/app/**

The default motif searched in this script is:
**Y..[YFT].{0,6}$**
