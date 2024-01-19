from Bio import SeqIO
import re
import numpy as np

# Function to extract gene name from the sequence header
def extract_gene_name(seq_record):
    match = re.search(r'gene:([^ ]+)', seq_record.description)
    return match.group(1)


# Function to extract genomic coordinates from the sequence header
def extract_coordinates(seq_record):
    # Check if "scaffold" is present in the header
    scaffold_match = re.search(r'supercontig:(\S+):(\S+):(\d+):(\d+):(-?\d+)', seq_record.description)
    if scaffold_match:
        _, scaffold, start, end, strand = scaffold_match.groups()
        return scaffold, int(start), int(end), int(strand)

    # Check if "chromosome" is present in the header
    chromosome_match = re.search(r'chromosome:(\S+):(\d+):(\d+):(\d+):(-?\d+)', seq_record.description)
    if chromosome_match:
        _, chromosome, start, end, strand = chromosome_match.groups()
        return chromosome, int(start), int(end), int(strand)
    return None

# Function to calculate the physical distance between two genes
def calculate_distance(gene_sequences, gene1_name, gene2_name):
    gene1_record = next((record for record in gene_sequences.values() if extract_gene_name(record) == gene1_name), None)
    coordinates_gene1 = extract_coordinates(gene1_record)
    scaffold_gene1, start_gene1, end_gene1, strand_gene1 = coordinates_gene1

    gene2_record = next((record for record in gene_sequences.values() if extract_gene_name(record) == gene2_name), None)
    coordinates_gene2 = extract_coordinates(gene2_record)
    scaffold_gene2, start_gene2, end_gene2, strand_gene2 = coordinates_gene2

    if scaffold_gene1 == scaffold_gene2:
        distance = abs(start_gene1 - start_gene2)
        return distance, scaffold_gene1, scaffold_gene2
    else:
        return np.nan, scaffold_gene1, scaffold_gene2

fasta_file = 'Medicago_truncatula.MedtrA17_4.0.cds.all.fa'
gene_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

input_file_path = 'dS_values_final'
output_file_path = 'dS_values_distances.txt'

with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
    for line in infile:
        gene1_name, gene2_name, ds = line.split()
        result = calculate_distance(gene_sequences, gene1_name, gene2_name)
        outfile.write(f"{gene1_name}\t{gene2_name}\t{ds}\t{result[0]}\t{result[1]}\t{result[2]}\n")

