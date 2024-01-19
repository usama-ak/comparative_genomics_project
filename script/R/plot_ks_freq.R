library(data.table)
library(mclust)

# Read data
dS_values <- fread("dS_values_distances.txt")
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