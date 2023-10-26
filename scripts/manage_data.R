# Loading datasets
A <- read.csv("data/SetA_lipid_FLO2019Control.csv")
B <- read.csv("data/SetB_lipid_FLO2022_lowP.csv")


###### Filtering just the Peak intensities 

# Find column names that match the pattern "S1_Run" at the end
selected_columns_A <- grep("S1_Run\\d+", colnames(A), value = TRUE)
selected_columns_B <- grep("S1_Run\\d+", colnames(B), value = TRUE)
#selected_columns <- grep("ISTD", colnames(A), value = TRUE)


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


#write.csv(A_filtered,"results/A_filtered.csv", row.names = T)
#write.csv(B_filtered,"results/B_filtered.csv", row.names = T)

##### Filtering based on lipids
lipid_A <- read.csv("data/lipid_class_A.csv")
lipid_A$X.Scan. <- as.character(lipid_A$X.Scan.)
lipid_B <- read.csv("data/lipid_class_B.csv")
lipid_B$X.Scan. <- as.character(lipid_B$X.Scan.)

# Filtering only PC's
# Create a logical vector indicating rows where Compound_Name starts with "PC"
#pc_names <- grepl("^PC", lipid$Compound_Name)

# Subset the data frame to include only rows with PC names
#pc_data <- lipid[pc_names, ]

# View the first few rows of the subsetted data frame
#str(pc_data)

#unique_scans <- pc_data %>%
#  group_by(X.Scan.) %>%
#  summarize(Compound_Name = toString(unique(Compound_Name)))


# Unique ID's  
unique_scans_all_lipid_A <- lipid_A %>%
  dplyr::select(Compound_Name,X.Scan.)
unique_scans_all_lipid_B <- lipid_B %>%
  dplyr::select(Compound_Name,X.Scan.)


# Group and summarize unique_scans_all_lipid to concatenate Compound_Name values
# result <- unique_scans_all_lipid %>%
#   dplyr::group_by(X.Scan.) %>%
#   dplyr::summarize(Compound_Name = paste(Compound_Name, collapse = ","))
# 
# # Merge the summarized data with A_filtered based on X.Scan.
# A_filtered_lipid <- inner_join(result, A_filtered, by = "X.Scan.")
# B_filtered_lipid <- inner_join(result, B_filtered, by = "X.Scan.")
# 
# 
# length(unique(unique_scans_all_lipid$X.Scan.))
# length(unique(A_filtered$X.Scan.))
# length(unique(B_filtered$X.Scan.))
# 
# length(intersect(unique(A_filtered$X.Scan.),unique(B_filtered$X.Scan.)))
# length(intersect(A_filtered_lipid$X.Scan.,B_filtered_lipid$X.Scan.))
# intersect(A_filtered_lipid$Compound_Name,B_filtered_lipid$Compound_Name)

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

# Subset rows in A_filtered using row names from unique_scans
subset_A <- inner_join(A_filtered, unique_scans_all_lipid_A)
subset_A <- subset_A %>%
  dplyr::select(-X.Scan.) %>%
  dplyr::select(Compound_Name, everything())
  
subset_B <- inner_join(B_filtered, unique_scans_all_lipid_B)
subset_B <- subset_B %>%
  dplyr::select(-X.Scan.) %>%
  dplyr::select(Compound_Name, everything())


# Subsetting just the PC's
# Extract the "PC" values from the selected column names
PC_A <- grep("PC", subset_A$Compound_Name, value = TRUE)
PC_B <- grep("PC", subset_B$Compound_Name, value = TRUE)

# Subsetting the PCs
subset_A_PC <- as.data.frame(subset_A[subset_A$Compound_Name %in% PC_A, ])
subset_B_PC <- as.data.frame(subset_B[subset_B$Compound_Name %in% PC_B, ])


