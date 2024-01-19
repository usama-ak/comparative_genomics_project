README Comparative Genomics
Last Name : AKHTAR
First Name : USAMA

***********************
***********************

Participation in :

Bash script that calculated the length of the proteins and add the length in a new column in the file 'Medicago_truncatila_liste'.
(See git for the script : '1_prot_gene_length.sh' )

Output File : Medicago_truncatula_protein_length

***********************
***********************

Create python script to preprocessing and merging data.

Function Merging Files to merge the blast result file and the protein length file generated before.

'''
import os

# Function to read files and merge data
def merge_files(file1_path, file2_path, output_path, a, b):
    # Create a dictionary to store the second and third column values of file2
    file2_data = {}
    
    with open(file2_path, 'r') as file2:
        for line in file2:
            columns = line.strip().split()
            if len(columns) >= 3:
                key = columns[0]
                value2 = columns[1]
                value3 = columns[2]
                file2_data[key] = (value2, value3)

    with open(file1_path, 'r') as file1, open(output_path, 'w') as output_file:
        for line in file1:
            columns = line.strip().split()
            if len(columns) >= a:
                key = columns[b]
                if key in file2_data:
                    value2, value3 = file2_data[key]
                    output_line = f"{line.strip()}\t{value2}\t{value3}\n"
                    output_file.write(output_line)
                else:
                    output_file.write(f"{line.strip()} \t \n")
                    
'''

Function Filtering Hits to filter the blast hit and retrieve the hits with the maximum bitscore for the alignments in both directions

'''
def filter_hits2(file1_path, output_path):
    mydict = {}

    # Read file1 and store bitscore values in a dictionary for both directions
    with open(file1_path, 'r') as file1:
        for line in file1:
            columns = line.strip().split()
            key1 = columns[0]
            key2 = columns[1]
            bitscore = float(columns[11])  # Convert bitscore to float

            # Consider both directions: (key1, key2) and (key2, key1)
            mydict[(key1, key2)] = max(mydict.get((key1, key2), -float('inf')), bitscore)
            mydict[(key2, key1)] = max(mydict.get((key2, key1), -float('inf')), bitscore)

    # Read file1 again, filter data based on bitscore, and write the filtered data to the output file
    with open(file1_path, 'r') as file1, open(output_path, 'w') as output_file:
        for line in file1:
            columns = line.strip().split()
            key = (columns[0], columns[1])

            # Retrieve the maximum bitscore for both directions
            max_bitscore = max(mydict.get((key[0], key[1]), -float('inf')),
                              mydict.get((key[1], key[0]), -float('inf')))

            if float(columns[11]) == max_bitscore:  # Only write if it's the maximum bitscore
                output_file.write(f"{line.strip()} \t \n")
'''
(See git for the script : '2_column_comp.py' )

Output file : filtered_data.txt

***********************
***********************

Filtering based on identity and coverage : 

30% identity and 50% coverage:
awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=30 && ($8-$7+1)/$14 >=0.5 && ($10-$9+1)/$16 >= 0.5 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_30_50.txt;

30% identity and 80% coverage:
awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=30 && ($8-$7+1)/$14 >=0.8 && ($10-$9+1)/$16 >= 0.8 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_30_80.txt;


Clustering : MCL command lines: 

30% identity and 50% coverage:
mcl mcl_input_30_50.txt --abc -o mcl_output_30_50;

30% identity and 80% coverage:
mcl mcl_input_30_80.txt --abc -o mcl_output_30_80;

Number of families and max family size:

awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_30_50;
Total Families: 6018 Max Family Size: 605

awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_30_80;
Total Families: 5213 Max Family Size: 285

(See git for the script : '3_mcl_processing.sh' )

Output files : mcl_output_30_50, mcl_output_30_80

***********************
***********************

Plot the family size with their frequency to make the best choice

(See git for the script : '4_families_plots.py' )

***********************
***********************

Random sampling to have less members per family using R.


file_path <- "mcl_output_30_50"
data <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, quote = "")

# Convert the data frame to a list of vectors
data_vector <- lapply(seq_len(nrow(your_data)), function(i) unname(unlist(your_data[i, ])))
data_vector <- lapply(data_vector, function(x) x[!is.na(x) & x != ""])

family_sizes <- sapply(data_vector, length)

# Number of members per family
members_per_family <- 100 

# Function to sample members from each family
sample_members <- function(family, n) {
  if (length(family) <= n) {
    return(family)
  } else {
    return(sample(family, n, replace = FALSE))
  }
}

# Apply the sampling function to each family in data_vector
sampled_data <- lapply(data_vector, sample_members, n = members_per_family)


file_path <- "output_families_30_50.txt"
file_conn <- file(file_path, "w")

# Write each vector as a line in the file
for (i in 1:length(sampled_data)) {
  vec <- paste(sampled_data[[i]], collapse = "\t")
  writeLines(vec, file_conn)
}
close(file_conn)

(See git for the script : '5_genes_pairs.py' )

Output files : output_families_30_50.txt

***********************
***********************

Generate pairs of genes within each family generated by MCL.

(See git for the script : '5_genes_pairs.py' )

Output file : gene_pairs_output_30_50.txt

***********************
***********************

Create script to follow the pipeline for PAML.

Output file : dS_values

***********************
***********************

We filtered the Ks values < 3 = dS_values_final file

Plot the frequency of Ks values with R 

library(data.table)
library(mclust)

# Read data
dS_values <- fread("dS_values_final")
dS_values_filtered <- dS_values[dS_values$V3 < 3, ]
ks_values <- as.numeric(dS_values_filtered$V3)

# Fit a Gaussian Mixture Model (GMM) with two components
gmm_model <- Mclust(ks_values, G = 2)

# Plot the histogram
hist_obj <- hist(ks_values, breaks = 150, main = "Distribution of Ks Values", 
                 xlab = "Ks", ylab = "Gene pairs", col = "lightblue", border = "black",
                 xlim = c(0, 3))

# Add GMM components to the plot
lines(density(ks_values), col = "blue", lty = 2, lwd = 2)

# Plot the density curve with the adjusted bandwidth
lines(density(ks_values, bw = sqrt(as.numeric(variance1))), col = "red", lty = 2, lwd = 2)
lines(density(ks_values, bw = sqrt(as.numeric(variance2))), col = "green", lty = 2, lwd = 2)

# Highlight the means of the GMM components
abline(v = mean1, col = 'red', lwd = 3)
abline(v = mean2, col = 'green', lwd = 3)

(See git for the script : 'plot_ks_freq.R' )

***********************
***********************

Student Test to see if the peaks are significative : 

library(data.table)
library(mclust)

# Read data
dS_values <- fread("dS_values_final")
dS_values_filtered <- dS_values[dS_values$V3 < 3, ]
ks_values <- as.numeric(dS_values_filtered$V3)

# Fit a Gaussian Mixture Model (GMM) with two components
gmm_model <- Mclust(ks_values, G = 2)

# Perform a t-test
t_test_result <- t.test(ks_values ~ gmm_model$classification)
p_value <- t_test_result$p.value
alpha <- 0.05

# Print the p-value
print(p_value)

# Check for significance
if (p_value < alpha) {
  cat("The difference in means is statistically significant.\n")
} else {
  cat("The difference in means is not statistically significant.\n")
}

(See git for the script : 'student_test.R' )

***********************
***********************

Compute the distance between the genes in each pair if they are in the same scaffold or chromosome, using python script.

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

with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
    for line in infile:
        gene1_name, gene2_name, ds = line.split()
        result = calculate_distance(gene_sequences, gene1_name, gene2_name)
        outfile.write(f"{gene1_name}\t{gene2_name}\t{ds}\t{result[0]}\t{result[1]}\t{result[2]}\n")

(See git for the script : 'calcul_distances.py' )

Output file : dS_values_distances.txt

***********************
***********************

Python script to plot the distances : bar plot and box plot

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv('dS_values_distances.txt', sep='\t', header=None, names=['Gene1', 'Gene2', 'Distance', 'Scaffold1', 'Scaffold2'])

# Add a new column 'SameChromosome' indicating whether genes are on the same chromosome
data['SameChromosome'] = data['Distance'].notna()

# Count the occurrences of genes on the same chromosome
same_chromosome_counts = data['SameChromosome'].value_counts()

# Bar plot
same_chromosome_counts.plot(kind='bar', color=['red', 'green'])
plt.xticks([0, 1], ['Different Chromosomes', 'Same Chromosome'], rotation=0)
plt.xlabel('Chromosome Relationship')
plt.ylabel('Count')
plt.title('Count of Genes on Different and Same Chromosomes')
plt.show()

# Filter out rows where 'Distance' is NaN
computed_distances = data[data['Distance'].notna()]

# Boxplot of distances for each scaffold
plt.figure(figsize=(12, 6))
ax = sns.boxplot(x='Scaffold1', y='Distance', data=computed_distances)
ax.set_yscale('log')
plt.xlabel('Scaffold')
plt.ylabel('Distance')
plt.title('Distribution of Distances for Each Scaffold')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

(See git for the script : 'distance_plots.py' )

***********************
***********************

Check the chromosomes involved in the pairs of genes in each peaks to see if we have all 8 of them.

library(data.table)
library(mclust)

# Read data
dS_values <- fread("dS_values_distances.txt")
dS_values_filtered <- dS_values[dS_values$V3 < 3, ]
ks_values <- as.numeric(dS_values_filtered$V3)

# Fit a Gaussian Mixture Model (GMM) with two components
gmm_model <- Mclust(ks_values, G = 2)

# Extracting relevant columns
chromosomes_df <- dS_values_filtered[, c("V5", "V6")]

# Add cluster assignments to the data frame
chromosomes_df$Cluster <- gmm_model$classification

# Define a function to print the chromosomes for a cluster
print_chromosomes <- function(cluster_data) {
  unique_chromosomes <- unique(unlist(cluster_data[, c("V5", "V6")]))
  cat("Cluster contains the following chromosomes:", paste(unique_chromosomes, collapse = ", "), "\n")
}

# Apply the function to each cluster
for (cluster in unique(chromosomes_df$Cluster)) {
  cluster_data <- chromosomes_df[chromosomes_df$Cluster == cluster, ]
  print_chromosomes(cluster_data)
  }


***********************
***********************

Now we wanted to see all the possible chromosomes involved in duplicated genes and count their occurences

library(data.table)
library(mclust)
library(ggplot2)

# Read data
dS_values <- fread("dS_values_distances.txt")
dS_values_filtered <- dS_values[dS_values$V3 < 3, ]
ks_values <- as.numeric(dS_values_filtered$V3)

# Fit a Gaussian Mixture Model (GMM) with two components
gmm_model <- Mclust(ks_values, G = 2)


chromosomes_df <- chromosomes_df[!(grepl("^scaffold", chromosomes_df$V5) | 
                                     grepl("^scaffold", chromosomes_df$V6) | 
                                     chromosomes_df$V5 == chromosomes_df$V6), ]

chromosome_counts <- table(chromosomes_df$V5, chromosomes_df$V6)

# Convert the table to a data frame
counts_df <- as.data.frame(chromosome_counts)
colnames(counts_df) <- c("Chromosome1", "Chromosome2", "Count")

# Plot the bar chart
ggplot(counts_df, aes(x = interaction(Chromosome1, Chromosome2, lex.order = TRUE), y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Chromosome Pairing", y = "Count", title = "Chromosome Pairing Counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

(See git for the script : 'chr_combinations.R' )


***********************
***********************

Finally we wanted to see the relationship between the distances and the 2 peaks we found : distribution of distances per peak.

library(data.table)
library(mclust)
library(ggplot2)

# Read data
dS_values <- fread("dS_values_distances.txt")
dS_values_filtered <- dS_values[dS_values$V3 < 3, ]
ks_values <- as.numeric(dS_values_filtered$V3)

# Fit a Gaussian Mixture Model (GMM) with two components
gmm_model <- Mclust(ks_values, G = 2)

# Extract cluster assignments
cluster_assignments <- gmm_model$classification

# Extract distances for each cluster
cluster_distances <- lapply(unique(cluster_assignments), function(cluster_num) {
  cluster_data <- dS_values[!is.na(dS_values$V4) & cluster_assignments == cluster_num, "V4"]
  return(cluster_data)
})

# Convert cluster_distances to a data frame for boxplot
cluster_distances_df <- data.frame(Cluster = rep(1:length(cluster_distances), sapply(cluster_distances, length)),
                                   Distance = unlist(cluster_distances))

# Plot boxplots
ggplot(cluster_distances_df, aes(x = factor(Cluster), y = Distance, fill = factor(Cluster))) +
  geom_boxplot() +
  labs(x = "Peaks", y = "Distance", title = "Boxplots of Distances for Each Peaks") +
  theme_minimal()
