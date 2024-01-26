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
selected_columns_A <- grep("S1_Run\\d+", colnames(A), value = TRUE)
selected_columns_B <- grep("S1_Run\\d+", colnames(B), value = TRUE)

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
columns_to_keep_A <- which(zero_proportions_A < threshold)
columns_to_keep_B <- which(zero_proportions_B < threshold)

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


################################################################################
########## Averaging the values of similar compounds row-wise  #################
################################################################################


# Remove "_<number>" from the end of the strings in 'Column1'
subset_A$Compound_Name <- gsub("_\\d+$", "", subset_A$Compound_Name)
subset_B$Compound_Name <- gsub("_\\d+$", "", subset_B$Compound_Name)


# Loop through all columns except the first one ('Column1')
for (col in names(subset_A)[-1]) {
  subset_A[[col]] <- as.numeric(subset_A[[col]])
  subset_A[[col]][is.na(subset_A[[col]])] <- 0  # Replace NA with 0 for non-numeric values
}
for (col in names(subset_B)[-1]) {
  subset_B[[col]] <- as.numeric(subset_B[[col]])
  subset_B[[col]][is.na(subset_B[[col]])] <- 0  # Replace NA with 0 for non-numeric values
}

# Group by 'Column1' and calculate the mean for other columns
averaged_subset_A <- aggregate(. ~ Compound_Name, data = subset_A, FUN = mean)
averaged_subset_B <- aggregate(. ~ Compound_Name, data = subset_B, FUN = mean)


################################################################################
########## Normalizing by centering to median of log10 values  #################
################################################################################

# Select only the numeric columns (excluding the first column)
numeric_subset_A <- averaged_subset_A[, -1]
numeric_subset_B <- averaged_subset_B[, -1]

# Log10 normalize the numeric columns
log10_normalized_subset_A <- log10(numeric_subset_A)
log10_normalized_subset_B <- log10(numeric_subset_B)

# Center by median for log10 normalized values
centered_subset_A <- apply(log10_normalized_subset_A, 2, function(x) x - median(x, na.rm = TRUE))
centered_subset_B <- apply(log10_normalized_subset_B, 2, function(x) x - median(x, na.rm = TRUE))

# Combine the first column with the centered log10 normalized values
normalized_subset_A <- data.frame(averaged_subset_A[, 1], centered_subset_A)
normalized_subset_B <- data.frame(averaged_subset_B[, 1], centered_subset_B)

# You can assign column names if needed
colnames(normalized_subset_A) <- c("Compound_Name", colnames(centered_subset_A))
colnames(normalized_subset_B) <- c("Compound_Name", colnames(centered_subset_B))

# Remove the unwanted variables
rm(list = setdiff(ls(), c("normalized_subset_A", "normalized_subset_B")))
