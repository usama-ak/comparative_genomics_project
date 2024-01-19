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
