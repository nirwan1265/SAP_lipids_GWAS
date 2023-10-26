# Find common compounds between PC_A and PC_B
common_compounds <- intersect(PC_A$Compound_Prefix, PC_B$Compound_Prefix)

# Subset the data frames to include only common compounds
common_compounds_A <- PC_A[PC_A$Compound_Prefix %in% common_compounds, ]
common_compounds_B <- PC_B[PC_B$Compound_Prefix %in% common_compounds, ]


str(common_compounds_A)


# Initialize p-values vector
p_values <- numeric(length(common_compounds))

# Initialize p-values vector
p_values <- numeric(nrow(common_compounds_A))

# Perform statistical tests for each compound row-wise
for (i in 1:nrow(common_compounds_A)) {
  data_A <- as.numeric(common_compounds_A[i, -1])
  data_B <- as.numeric(common_compounds_B[i, -1])
  
  # Perform a statistical test (e.g., Wilcoxon rank-sum test) for this compound
  test_result <- wilcox.test(data_A, data_B)
  p_values[i] <- test_result$p.value
}

# Add the Compound_Prefix column to the p-values vector
p_values <- cbind(Compound_Prefix = common_compounds_A$Compound_Prefix, P_Value = p_values)

str(p_values)
# Convert p_values into a data frame
p_values <- as.data.frame(p_values, stringsAsFactors = FALSE)

# Sort the p-values data frame by P_Value (ascending order)
p_values <- p_values[order(as.numeric(p_values$P_Value)), ]


# Differences between the concentration
# Create a new dataframe to store the differences
differences <- common_compounds_A

# Exclude the 'Compound_Prefix' column from both dataframes
common_compounds_A_values <- common_compounds_A[, -1]
common_compounds_B_values <- common_compounds_B[, -1]

# Calculate row-wise averages for common_compounds_A and common_compounds_B
average_A <- rowMeans(common_compounds_A[, -1])
average_B <- rowMeans(common_compounds_B[, -1])

# Create a new dataframe to store the differences
differences <- data.frame(Compound_Prefix = common_compounds_A$Compound_Prefix, Difference = average_A - average_B)

# Find unique compounds in data frame A
unique_compounds_A <- setdiff(PC_A$Compound_Prefix, PC_B$Compound_Prefix)

# Find unique compounds in data frame B
unique_compounds_B <- setdiff(PC_B$Compound_Prefix, PC_A$Compound_Prefix)
