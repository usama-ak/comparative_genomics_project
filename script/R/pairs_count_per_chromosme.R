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

