################################################################################
############ Loading the Quantitative Peaks Intensities.  ######################
################################################################################

# Control
A <- read.csv("data/SetA_lipid_FLO2019Control.csv")

# Low P
B <- read.csv("data/SetB_lipid_FLO2022_lowP.csv")


################################################################################
############## Filtering just the Sample Peak intensities ######################
################################################################################

# Find column names that match the pattern "S1_Run" at the end
selected_columns_A <- grep("ISTD", colnames(A), value = TRUE)
selected_columns_B <- grep("ISTD", colnames(B), value = TRUE)

# Select the corresponding columns in your data frame A
A_filtered <- A[, colnames(A) %in% selected_columns_A]
B_filtered <- B[, colnames(B) %in% selected_columns_B]

# Extract the "PI" values from the selected column names
pi_values_A <- paste0("PI",sub(".*_PI(\\d+)_.*", "\\1", selected_columns_A))
pi_values_B <- paste0("PI",sub(".*_PI(\\d+)_.*", "\\1", selected_columns_B))

# Rename the columns in A_filtered with the cleaned "PI" values
colnames(A_filtered) <- pi_values_A
colnames(B_filtered) <- pi_values_B

# Adding colnames
rownames(A_filtered) <- A[,1]
rownames(B_filtered) <- B[,1]
rm(A,B)

# Removing CHECKs
A_filtered <- A_filtered[, !grepl("CHECK", names(A_filtered))]
B_filtered <- B_filtered[, !grepl("Check", names(B_filtered))]


################################################################################
############ Filtering out samples with more than 50% zeroes ###################
################################################################################

# Set the threshold for proportion of zeroes you want to remove
threshold <- 0.5

# Calculate the proportion of zeroes in each column
zero_proportions_A <- colMeans(A_filtered == 0, na.rm = TRUE)
zero_proportions_B <- colMeans(B_filtered == 0, na.rm = TRUE)

# Get the column indices to keep (where less than 50% are zeroes)
columns_to_keep_A <- which(zero_proportions_A > threshold)
columns_to_keep_B <- which(zero_proportions_B > threshold)

# Subset A_filtered to include only the selected columns
A_filtered <- A_filtered[, columns_to_keep_A]
A_filtered$X.Scan. <- rownames(A_filtered)
B_filtered <- B_filtered[, columns_to_keep_B]
B_filtered$X.Scan. <- rownames(B_filtered)


################################################################################
###################### Loading the Lipid Datasets  #############################
################################################################################


# Control
lipid_A <- read.csv("data/lipid_class_A.csv")
lipid_A$X.Scan. <- as.character(lipid_A$X.Scan.)

# Low P
lipid_B <- read.csv("data/lipid_class_B.csv")
lipid_B$X.Scan. <- as.character(lipid_B$X.Scan.)

# Selecting only the requried columns
unique_scans_all_lipid_A <- lipid_A %>%
  dplyr::select(Compound_Name,X.Scan.)
unique_scans_all_lipid_B <- lipid_B %>%
  dplyr::select(Compound_Name,X.Scan.)

rm(lipid_A,lipid_B)

# Process duplicates
unique_scans_all_lipid_A$Compound_Name <- with(unique_scans_all_lipid_A, ave(Compound_Name, Compound_Name, FUN = function(x) {
  if (length(x) > 1) {
    return(paste0(x, "_", 1:length(x)))
  } else {
    return(x)
  }
}))

unique_scans_all_lipid_B$Compound_Name <- with(unique_scans_all_lipid_B, ave(Compound_Name, Compound_Name, FUN = function(x) {
  if (length(x) > 1) {
    return(paste0(x, "_", 1:length(x)))
  } else {
    return(x)
  }
}))

# Convert X.Scan. column to character type
unique_scans_all_lipid_A$X.Scan. <- as.character(unique_scans_all_lipid_A$X.Scan.)
unique_scans_all_lipid_B$X.Scan. <- as.character(unique_scans_all_lipid_B$X.Scan.)


################################################################################
####################### Getting the Lipid Names  ###############################
################################################################################

# Subset rows in A_filtered using row names from unique_scans
subset_A <- inner_join(A_filtered, unique_scans_all_lipid_A)
subset_A <- subset_A %>%
  dplyr::select(-X.Scan.) %>%
  dplyr::select(Compound_Name, everything())

subset_B <- inner_join(B_filtered, unique_scans_all_lipid_B)
subset_B <- subset_B %>%
  dplyr::select(-X.Scan.) %>%
  dplyr::select(Compound_Name, everything())

rm(unique_scans_all_lipid_A,unique_scans_all_lipid_B)

################################################################################
####### Changing 0's in row to 2/3 lowest non-zero value row-wise  #############
################################################################################

# Function to replace zeros with 2/3 of the lowest non-zero value
replace_zeros <- function(row) {
  numeric_values <- as.numeric(row[-1])  # Convert numeric columns to numeric
  non_zero_values <- numeric_values[numeric_values != 0]  # Extract non-zero values
  if (length(non_zero_values) > 0) {
    lowest_non_zero <- min(non_zero_values)  # Find the lowest non-zero value
    row[-1][numeric_values == 0] <- 2/3 * lowest_non_zero  # Replace zeros with 2/3 of the lowest non-zero value
  }
  return(row)
}

# Identify numeric columns (excluding the first non-numeric column)
numeric_cols_A <- sapply(subset_A[-1], is.numeric)
numeric_cols_B <- sapply(subset_B[-1], is.numeric)

# Apply the replace_zeros function to numeric columns only
subset_A[, numeric_cols_A] <- t(apply(subset_A[, numeric_cols_A], 1, replace_zeros))
subset_B[, numeric_cols_B] <- t(apply(subset_B[, numeric_cols_B], 1, replace_zeros))



