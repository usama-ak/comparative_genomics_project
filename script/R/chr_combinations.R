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
