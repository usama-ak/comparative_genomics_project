library(data.table)
library(mclust)

# Read data
dS_values <- fread("dS_values_final")
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

