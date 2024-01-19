#!/usr/bin/python3

def generate_gene_pairs(input_file, output_file):
    gene_families = []
    
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Process each line as one family
    for line in lines:
        genes = line.strip().split('\t')
        gene_families.append(genes)

    # Generate gene pairs within families
    gene_pairs = []
    for family in gene_families:
        for i in range(len(family)):
            for j in range(i + 1, len(family)):
                gene_pairs.append((family[i], family[j]))

    # Write gene pairs to the output file
    with open(output_file, 'w') as output_file:
        for pair in gene_pairs:
            output_file.write('\t'.join(pair) + '\n')


input_file_path = 'output_families_30_50.txt'
output_file_path = 'gene_pairs_output_30_50.txt'

generate_gene_pairs(input_file_path, output_file_path)


