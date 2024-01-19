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










