README module Comparative Genomics 2023-2024 Project part

LAST NAME : YOUNCHA
FIRST NAME : Ali

***************************
*  Friday 20 october 2023 *
***************************

Work of the day :

The first step of our analysis was to retreive, the different proteins length (that are found in the BLAST output file) and then to add the corresponding value to each row of the file containing the list of protein IDs and their corresponding gene IDs using the following bash script :

#!/usr/bin/bash

## This bash script is able to find the protein identifier in the Medicago_truncatula.MedtrA17_4.0.pep.all.fa file, retrieves its sequence and calculates it's length.

awk '/^>/ { 
    if (length(protein) > 0) { 
        print length(protein); 
    }
    protein = ""; 
    next; 
}
{ 
    protein = protein $0; 
}
END { 
    if (length(protein) > 0) { 
        print length(protein); 
    }
}' Medicago_truncatula.MedtrA17_4.0.pep.all.fa > protein_length;
sed '1d' Medicago_truncatula_liste > new_list
paste new_list protein_length > Medicago_truncatula_protein_length;
rm protein_length;
rm new_list;

Then, I participated in developing a python script that was able to read the Medicago truncatula BLAST file results, and the file that was previously generated with the protein_ID, Gene_ID and the length of the protein. The idea here will be to add, to the BLAST file of all the proteins of the species Medicago Truncatula, the Gene accession ID and the corresponding protein length. We were able to add 4 more columns on the file (Gene_ID, protein length of the query, Gene_ID and protein length of the subject, respectively). Then it will filter hits (For example if we have A >> B | B >> A) to keep only the alignment that has the highest bit score so that we obtain reliable results at the end.

(See git for the script : '2_column_comp.py')

***************************
*  Friday 27 october 2023 *
***************************

Work of the day :

I performed a second filtering on the resulting file of the previous step by taking into account different percents of identity and coverage, then I performed MCL (The Markov Cluster Algorithm) on different datasets to generate genes that are clustered together (and therefore the corresponding family) knowing that each cluster represents a family. I only did these steps on the last two datasets.

To calculate the percent of coverage, I used the following awk commands :

awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=50 && ($8-$7+1)/$14 >=0.5 && ($10-$9+1)/$16 >= 0.5 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_50_50.txt;
awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=50 && ($8-$7+1)/$14 >=0.8 && ($10-$9+1)/$16 >= 0.8 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_50_80.txt;


Then, the following commands to run MCL :

mcl mcl_input_50_50.txt --abc -o mcl_output_50_50;
mcl mcl_input_50_80.txt --abc -o mcl_output_50_80;

And the following commands to retrieve the total number :


awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_50_50;
awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_50_80;


Then, I developed a python script that will plot, for each dataset the length of families on the X abciss and the frequency (number of families) on the Y abciss for all the dataset using the following script :

## This script will go through all the files and plot the length of the families based on the number of families that were observed for each length.

#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np

file_paths = ['mcl_output_30_50', 'mcl_output_30_80', 'mcl_output_50_50', 'mcl_output_50_80']

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
axes = axes.flatten()

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

(see Figure 10 on the report file)

***************************
*  Monday 25 december 2023 *
***************************

Work of the day :

Today, I developed a python script that will go through the output file of the MCL method, and will generate all the possible pairs of genes for each family.

 
#!/usr/bin/python3

def generate_gene_pairs(input_file, output_file):

    gene_families = []
    

    with open(input_file, 'r') as file:
        lines = file.readlines()


    for line in lines:
        genes = line.strip().split('\t')
        gene_families.append(genes)

    gene_pairs = []


    for family in gene_families:
        # Iterating gene_families variable using a for loop
        for i in range(len(family)):
            for j in range(i + 1, len(family)):
                # Creating tuples representing a gene pair
                gene_pairs.append((family[i], family[j]))


    with open(output_file, 'w') as output_file:
        for pair in gene_pairs:
            output_file.write('\t'.join(pair) + '\n')


input_file_path = 'output_families_30_50.txt'
output_file_path = 'gene_pairs_output_30_50'


generate_gene_pairs(input_file_path, output_file_path)


***************************
*  Saturday 30 december 2023 *
***************************


Work of the day :


I participated in performing sampling on our file containing the possible list of pairs of genes, since that I obtained more than 700.000 pairs of genes, the idea was to keep information only for the families having a total size of 100 or less.

(see git for the script : sampling.R)

Then once again, I executed the code that generates all the possible pairs of genes, as a result I obtained a total number of 364376 pairs.

***************************
*  Friday 13 January 2024 *
***************************


Work of the day :


I participated in developing a python script that will run the PAML analysis on all the possible pairs of genes 2 by 2. The objective of this step is to calculate and extract the Ks value of each of the pairs of the genes, to have at the end a file containing the gene pairs (column 1 and 2) and the corresponding Ks value (column 3).

(see git for the script : 6_paml_process.py)

Since that the number of pairs of genes was very high, we devided the execution of the script into 5 laptops, this took a total of ~40h to calculate all the Ks values.

***************************
*  Monday 15 January 2024 *
***************************

Work of the day :


By the time the PAML process finished, I pursued the analysis by filtering the output of the PAML program with a value of Ks < 3, the resulting file contained only the pairs of genes that have a Ks value less than 3. Then I plotted the distribution of Ks values based on  density. I also did the plots for a value of Ks < 5

library(data.table)
library(mclust)

# Read data
dS_values <- fread("dS_values_distances.txt")
dS_values_filtered <- dS_values[dS_values$V3 < 3, ]
ks_values <- as.numeric(dS_values_filtered$V3)

# Fit a Gaussian Mixture Model (GMM) with two components
gmm_model <- Mclust(ks_values, G = 2)

# Plot the density curve
plot(density(ks_values), main = "Distribution of Ks Values", 
     xlab = "Ks", ylab = "Density", col = "blue", lty = 2, lwd = 2)

# Highlight the means of the GMM components
abline(v = gmm_model$parameters$mean[1], col = 'red', lwd = 3)
abline(v = gmm_model$parameters$mean[2], col = 'green', lwd = 3)


***************************
*  Thursday 18 January 2024 *
***************************

Work of the day :


I participated in the realization of the student t-test to check if the highest peaks that were obtained in the Ks distribution are significative or not.

(see git for the script : student_test.R)


Then, I also participated in calculating the distances for the pairs of the genes that are located in the same chromosome, also those who are located in different chromosomes. 

(see git for the script : calcul_distances.py)

For the pairs taking part of the same chromosomes, we were able to calculate the distance between the genes by taking into account the start of the second gene minus the start of the first gene. For those who are in different chromosomes, the results were represented as Nan.

Moreover, I participated then in developing a script to plot :

	- The total count of pairs of genes that are found in different chromosomes and the same chromosomes
	- The distribution of the distances among all the scaffolds

(see git for the script : distance_plots.py)
(see plots in the report file)

Then I developed a script in R to plot, the count of pairs of genes that were located in each chromosome (in other words, the pairs of genes that take part of each chromosome). I observed that chromosome 3 has the highest count of pairs of genes (59000) followed by chromosome 4 (52000) pairs.

library(ggplot2)

# Rename columns for better clarity
colnames(dS_values) <- c("Gene1", "Gene2", "kS", "distance", "Chromosome1", "Chromosome2")

# Combine both chromosome columns into a single vector
all_chromosomes <- c(dS_values$Chromosome1, dS_values$Chromosome2)

# Filter out scaffolds
chromosomes <- all_chromosomes[grepl("^\\d+$", all_chromosomes)]

# Create a data frame with the chromosome counts
chromosome_counts <- data.frame(Chromosome = names(table(chromosomes)), Count = as.numeric(table(chromosomes)))

# Plot using ggplot2
ggplot(chromosome_counts, aes(x = Chromosome, y = Count, fill = Chromosome)) +
  geom_bar(stat = "identity") +
  labs(title = "Gene Pairs Count by Chromosome",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()


Finally, I participated in developing once again, two different R scripts :

	- The first script plotting the chromosome pairing counts
	- The second script showing the boxplots of the distribution of distances between pairs of genes in both peaks that were previously obtained in the Ks distribution

(see git for the scripts : chr_combinations.R and distances_peaks.R)
(see Supp. Fig 5 and 7 on the report file)


***************************
*  Friday 19 January 2024 *
***************************
Work of the day :

Participated in the github creation to submit all the scripts and commands that were used during this project, and finished writing the report (Discussion and Conclusion parts).
