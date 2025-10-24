# Finding motifs in a FASTA file
'''Here we are given a fasta file containing amino acid sequences from Trypanosoma brucei 'Tgondii_Proteins.fasta' and txt file 'Tgondii_product_descriptions.txt
containing gene id and product description.

We are to find all occurrences of a given motif in these sequences and then obtain the 
        gene id , product description and the last 30 amino acid sequences and write it to another output file '''
        
print("Program to find motifs in a FASTA file and extract relevant information.")

# import necessary libraries
import re
import sys

# Define file paths
if len(sys.argv) != 4:
   print('Usage : python finding_motif_fasta_txt.py <fasta_file> <description_file> <output_file>')
   sys.exit(1)

fasta_file = sys.argv[1]
des_file = sys.argv[2]
output_file = sys.argv[3]

motif = r'Y..[YFT].{0,6}$' # Define the motif pattern

# Read product descriptions into a dictionary
product_dict = {}
try:
    with open(des_file, 'r') as df:
        next(df)  # Skip header line
        for line in df:
            parts = line.strip().split('\t')
            gene_id = parts[0]
            description = parts[2]
            product_dict[gene_id] = description
except Exception as e:
    print(f"Error reading {des_file}: {e}")
    sys.exit(1)


# Read the FASTA file 
# Write the information to the output file
try:
    with open(output_file, 'w') as of:
        of.write("Gene ID\tProduct Description\tLast 30 Amino Acids\n")
        with open(fasta_file) as ff:
            for line in ff:
                if line.startswith('>'):
                    gene = line.rstrip().lstrip('>')
                    seq = next(ff).rstrip()
                match = re.search(motif, seq)
                if match:
                    seq_id = gene
                    out_seq = seq[len(seq)-30:]
                    assert len(out_seq) == 30
                    match_seq = match.group()
                    out_seq_lower = out_seq.replace(match_seq, match_seq.lower())
                    description = product_dict.get(seq_id, "Unknown")
                    of.write(f"{seq_id}\t{description}\t{out_seq_lower}\n")
except Exception as e:

    print(f"Erro processing files: {e}")

    sys.exit(1)

print(f"Motif search complete. Results written to {output_file}.")