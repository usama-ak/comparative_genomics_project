library(data.table)
library(mclust)

# Read data
dS_values <- fread("dS_values_distances.txt")
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