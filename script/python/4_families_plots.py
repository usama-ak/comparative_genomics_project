#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np

file_paths = ['mcl_output_30_50', 'mcl_output_30_80', 'mcl_output_50_50', 'mcl_output_50_80']

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
axes = axes.flatten()

# Read data from each file and calculate total genes
total_genes = 0
for i, (file_path, ax) in enumerate(zip(file_paths, axes), start=1):
    with open(file_path, 'r') as file:
        data = [line.strip().split('\t') for line in file]
        family_sizes = [len(line) for line in data]

        # Calculate total genes for each file
        total_genes_in_file = sum(family_sizes)
        print(f'Total genes in File {file_path}: {total_genes_in_file}')
        total_genes += total_genes_in_file

        # Plot histogram for each file
        ax.hist(family_sizes, bins=np.arange(1, max(family_sizes) + 2), alpha=0.7)
        ax.set_title(f'{file_path}')
        ax.set_xlabel('Family Size')
        ax.set_ylabel('Frequency')

plt.tight_layout()
plt.show()




