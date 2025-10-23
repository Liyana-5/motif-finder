# problem 1
import re
import sys

# uncomment these line to try the command line
#if len(sys.argv) != 4:
#    print('Command line not complete. Please check and try again\n')

#input_fasta, input_table, output_table = sys.argv[1:]

# comment these lines to try the command line
input_table = '/home/kathryn/BIOL4292/Lab5/Tgondii_product_descriptions.txt'
input_fasta = '/home/kathryn/BIOL4292/Lab5/Tgondii_Proteins.fasta'
output_table = '/home/kathryn/BIOL4292/Lab5/Tgondii_endosomal_proteins.txt'

# the regex for the YXXphi motif
# look at the lecture recording if you don't understand this
pattern = r'Y..[YFT].{0,6}$'

# testing some cases where I should see a match
# all three options for FTY and extremes of the range
assert re.search(pattern, 'AAAAAAYRVT')
assert re.search(pattern, 'AAAAAAYRVF')
assert re.search(pattern, 'AAAAAAYRVY')
assert re.search(pattern, 'AAAAAAYRVTAAAAAA')

# test some cases where I should not see a match
# match is in sequence but not in the right place
assert not re.search(pattern, 'AAAAAAYRVTAAAAAAAA')
# close but inexact match
assert not re.search(pattern, 'AAAAAAYRVVAAAAAA')

# make a product descriptions dict before I start
product_descs = {}
try:
    with open(input_table) as input_table:
        # skip the header row - we don't need it
        next(input_table)
        for line in input_table:
            line = line.rstrip().split('\t')
            # I don't need the identifier in the second column
            gene = line[0]
            desc = line[2]
            product_descs[gene] = desc
# if we can't find this file, we're a bit stuck - exit!
except FileNotFoundError:
    sys.exit(f'File {input_table} not found, please try again\n')


# now read the fasta file
# we'll write the output in this loop, so open this first
with open(output_table, 'w') as output:
    try:
        with open(input_fasta) as fasta:
            for line in fasta:
                # defline should start with >
                if line.startswith('>'):
                    gene = line.rstrip().lstrip('>')
                    # if the previous line was a defline, the next line should be sequence
                    seq = next(fasta).rstrip()
                    # test if the sequence matches our regex
                    match = re.search(pattern, seq)
                    if match:

                        # if it does, we want to display only the last 30 aa in the output
                        output_seq = seq[len(seq) - 30 : ]
                        # check I got my indexing above right!
                        assert len(output_seq) == 30

                        # get the sequence of the match
                        # in the protein sequence, replace the match with the match in lower case
                        match_seq = match.group()
                        output_seq = output_seq.replace(match_seq, match_seq.lower())

                        # if the gene has a product description, retrieve it and write to file
                        if gene in product_descs:
                            desc = product_descs[gene]
                            output.write(f'{gene}\t{desc}\t{output_seq}\n')
    # if we can't find the input file, we're a bit stuck. Exit!
    except FileNotFoundError:
        sys.exit(f'File {input_fasta} not found, please try again\n')




# problem 2
# usually I would put this at the top of the file
# here, I have two scripts in one file, so I've put it here for clarity
import sys

# I haven't put a command line in this script because it's too confusing when there are two
# scripts in the same file.
# You could create one here, the same way as above.


# function to calculate rpk
# This allows us to test this function
# We can also do error handling in the function
# Note I have two except clauses to handle different types of error
def calc_rpk(count, length):
    try:
        rpk = count / (length / 1000)
    except ZeroDivisionError:
        rpk = 0
    except TypeError:
        sys.exit(f'Count and Length must be numerical values. We saw count {count} and length {length}.\n')
    return rpk

# some simple tests of the function
# the second tests the case when the denominator is 0
assert calc_rpk(1000, 100) == 10000
assert calc_rpk(1000, 0) == 0

# I will need the gene lengths for the rpk calculation
# get these first and make a dict
gene_lengths = {}
with open('/home/kathryn/BIOL4292/Lab5/Trypanosoma_brucei.fasta') as input_fasta:
    for line in input_fasta:
        # parse the fasta file as usual
        if line.startswith('>'):
            gene = line.rstrip().lstrip('>')
            seq = next(input_fasta).rstrip()
            # put the length of the sequence into the dict - we don't need the seq itself
            gene_lengths[gene] = len(seq)

# variables to populate
rpk_sum = 0 # sum of rpk values
gene_rpks = {} # dict of gene: rpk
with open('/home/kathryn/BIOL4292/Lab5/Trypanosoma_brucei_counts.txt') as counts_input:
    for line in counts_input:
        # two-column tab delimited, so split at tab and assign to variables
        gene, count = line.rstrip().split('\t')
        # for each gene, check that we have the length
        if gene in gene_lengths:
            # if we do, retrieve it and calculate rpk
            length = gene_lengths[gene]
            rpk = int(count) / (length / 1000)

            # add the rpk to the sum, and the gene: rpk pair to the dict
            rpk_sum += rpk
            gene_rpks[gene] = rpk

# tpm_sum will be used for testing
tpm_sum = 0
# we can't calculcate TPM until we have the final rpk_sum
# so, we must do this in a separate loop
# this is a case where we do need to iterate the dict
with open('/home/kathryn/BIOL4292/Lab5/Trypanosoma_brucei_TPM.txt', 'w') as output:
    # here we iterate the dict
    for gene, rpk in gene_rpks.items():
        # calc tpm and write to file
        tpm = rpk / (rpk_sum / 1000000)
        output.write(f'{gene}\t{tpm}\n')
        tpm_sum += tpm

# the sum should equal 1 million
# there is a rounding error here 
# I can still use assert to check that it's close (and to define how close I want to be)
assert 1000000 - tpm_sum < 0.0001

# if you use the unittest module we can user assertAlmostEqual to do this: https://docs.python.org/3/library/unittest.html