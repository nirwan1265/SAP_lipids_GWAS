# Bootstrap test
#https://stats.stackexchange.com/questions/25738/is-a-log-transformation-a-valid-technique-for-t-testing-non-normal-data

normalized_subset_A <- as.data.frame(t(normalized_subset_A))
colnames(normalized_subset_A) <- normalized_subset_A[1,]
normalized_subset_A <- normalized_subset_A[-1,]

normalized_subset_B <- as.data.frame(t(normalized_subset_B))
colnames(normalized_subset_B) <- normalized_subset_B[1,]
normalized_subset_B <- normalized_subset_B[-1,]



# Find common compounds between PC_A and PC_B
common_compounds <- intersect(colnames(normalized_subset_A), colnames(normalized_subset_B))

# Subset the data frames to include only common compounds
common_compounds_A <- normalized_subset_A[,common_compounds]
common_compounds_B <- normalized_subset_B[,common_compounds]

# Loop through all columns except the first one ('Column1')
for (col in names(common_compounds_A)) {
  common_compounds_A[[col]] <- as.numeric(common_compounds_A[[col]])
  common_compounds_A[[col]][is.na(common_compounds_A[[col]])] <- 0  # Replace NA with 0 for non-numeric values
}
for (col in names(common_compounds_B)) {
  common_compounds_B[[col]] <- as.numeric(common_compounds_B[[col]])
  common_compounds_B[[col]][is.na(common_compounds_B[[col]])] <- 0  # Replace NA with 0 for non-numeric values
}

(common_compounds_B)



# Step 2: Find unique columns in each data frame
unique_columns_A <- setdiff(colnames(common_compounds_A), colnames(common_compounds_B))
unique_columns_B <- setdiff(colnames(common_compounds_B), colnames(common_compounds_A))


install.packages("MKinfer")
library(MKinfer)

# Initialize an empty list to store the test results
bootstrap_results <- list()

# Loop through each column in the data frames
for (col in colnames(common_compounds_A)) {
  # Extract the columns as vectors
  x <- common_compounds_A[[col]]
  y <- common_compounds_B[[col]]
  
  # Check if both columns are numeric
  if (is.numeric(x) && is.numeric(y)) {
    # Perform the bootstrap t-test and store the results
    bootstrap_results[[col]] <- boot.t.test(x, y,
                                            alternative = "two.sided", # Change if needed
                                            mu = 0, # Hypothesized difference in means
                                            paired = FALSE, # Should be FALSE for unpaired test
                                            var.equal = FALSE, # Change to TRUE if variances are assumed equal
                                            conf.level = 0.95, # Confidence level
                                            R = 9999, # Number of bootstrap replicates
                                            symmetric = FALSE) # Change to TRUE to assume symmetry
  } else {
    # If the columns are not numeric, store a message instead of a result
    bootstrap_results[[col]] <- "Column is not numeric"
  }
}

# View the results
str(bootstrap_results[1])


# Initialize an empty dataframe to store the results
results <- data.frame(Compound = character(), PValue = numeric(), stringsAsFactors = FALSE)

# Loop through the bootstrap_results list
for (name in names(bootstrap_results)) {
  # Extract the p-value for the current list item
  p_value <- bootstrap_results[[name]]$p.value
  
  # Append to the results dataframe
  results <- rbind(results, data.frame(Compound = name, PValue = p_value))
}

# View the results dataframe
results

# Assuming you have the 'results' dataframe from the previous step
# Define the significance level, commonly 0.05
significance_level <- 0.05

# Filter the results to include only statistically significant p-values
sig_results <- results[results$PValue <= significance_level, ]

# View the statistically significant results dataframe
sig_results

# Subset the data frames to include only columns with significant p-values
sig_common_compounds_A <- common_compounds_A[, sig_results$Compound]
sig_common_compounds_B <- common_compounds_B[, sig_results$Compound]

# Check the structure of the subsetted data frames
str(sig_common_compounds_A)
str(sig_common_compounds_B)


# Find the common columns that are significant and present in both data frames
common_sig_columns <- intersect(colnames(common_compounds_A), colnames(common_compounds_B))
common_sig_columns <- common_sig_columns[common_sig_columns %in% sig_results$Compound]

# Subset the data frames to include only the common significant columns
sig_common_compounds_A <- common_compounds_A[, common_sig_columns, drop = FALSE]
sig_common_compounds_B <- common_compounds_B[, common_sig_columns, drop = FALSE]

# Check the structure of the subsetted data frames
str(sig_common_compounds_A)
str(sig_common_compounds_B)


# 
# common_compounds <- intersect(colnames(PC_A), colnames(PC_B))
# 
# # Subset the data frames to include only common compounds
# common_compounds_A <- PC_A[,common_compounds]
# common_compounds_B <- PC_B[,common_compounds]
# 
# 
# ### Log10 transform and centered around the median. 
# # Log10 transform the data frames common_compounds_A and common_compounds_B
# common_compounds_A_log10 <- log10(common_compounds_A)
# common_compounds_B_log10 <- log10(common_compounds_B)
#
# # Calculate the median for each column in the log10-transformed data frames
# median_A <- apply(common_compounds_A_log10, 2, median)
# median_B <- apply(common_compounds_B_log10, 2, median)
# 
# # Center the log10-transformed data frames around the median
# common_compounds_A <- common_compounds_A_log10 - median_A
# common_compounds_B <- common_compounds_B_log10 - median_B
# 
# 
# str(common_compounds_A)
# str(common_compounds_B)
# 
# 
# #Check if column names match
# if(all(names(common_compounds_A) == names(common_compounds_B))) {
#   # Initialize a data frame to store test results
#   results <- data.frame(Column = character(0), p_value = numeric(0))
#   
#   # Loop through column names and perform t-tests
#   for (col_name in names(common_compounds_A)) {
#     t_test_result <- t.test(common_compounds_A[[col_name]], common_compounds_B[[col_name]])
#     results <- rbind(results, data.frame(Column = col_name, p_value = t_test_result$p.value))
#   }
#   
#   # Print or analyze the results
#   print(results)
# } else {
#   print("Column names do not match between the two data frames.")
# }
# 
# # Sort the results data frame by p-value in ascending order
# sorted_results <- results[order(results$p_value), ]
# 
# # View the sorted results
# print(sorted_results)
# 
# # Significance level (adjust as needed)
# significance_level <- 0.05
# 
# # Filter the sorted results to include only statistically significant columns
# significant_results <- sorted_results[sorted_results$p_value < significance_level, ]
# 
# # View the significant results
# print(significant_results)
# 
# # Loop through the significant columns and calculate the mean difference
# for (i in 1:nrow(significant_results)) {
#   col_name <- significant_results$Column[i]
#   mean_diff <- mean(common_compounds_A[[col_name]]) - mean(common_compounds_B[[col_name]])
#   significant_results$Mean_Difference[i] <- mean_diff
# }
# 
# # View the updated significant_results data frame
# print(significant_results)
