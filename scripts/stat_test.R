# Find common compounds between PC_A and PC_B
common_compounds <- intersect(colnames(PC_A), colnames(PC_B))

# Subset the data frames to include only common compounds
common_compounds_A <- PC_A[,common_compounds]
common_compounds_B <- PC_B[,common_compounds]


str(common_compounds_A)
str(common_compounds_B)


#Check if column names match
if(all(names(common_compounds_A) == names(common_compounds_B))) {
  # Initialize a data frame to store test results
  results <- data.frame(Column = character(0), p_value = numeric(0))
  
  # Loop through column names and perform t-tests
  for (col_name in names(common_compounds_A)) {
    t_test_result <- t.test(common_compounds_A[[col_name]], common_compounds_B[[col_name]])
    results <- rbind(results, data.frame(Column = col_name, p_value = t_test_result$p.value))
  }
  
  # Print or analyze the results
  print(results)
} else {
  print("Column names do not match between the two data frames.")
}

# Sort the results data frame by p-value in ascending order
sorted_results <- results[order(results$p_value), ]

# View the sorted results
print(sorted_results)

# Significance level (adjust as needed)
significance_level <- 0.05

# Filter the sorted results to include only statistically significant columns
significant_results <- sorted_results[sorted_results$p_value < significance_level, ]

# View the significant results
print(significant_results)

# Loop through the significant columns and calculate the mean difference
for (i in 1:nrow(significant_results)) {
  col_name <- significant_results$Column[i]
  mean_diff <- mean(common_compounds_A[[col_name]]) - mean(common_compounds_B[[col_name]])
  significant_results$Mean_Difference[i] <- mean_diff
}

# View the updated significant_results data frame
print(significant_results)
