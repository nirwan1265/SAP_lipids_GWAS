# Loading data
A <- read.csv("data/SetA_lipid_FLO2020Control.csv")

###### Filtering just the Peak intensities 

# Find column names that match the pattern "S1_Run" at the end
selected_columns <- grep("S1_Run\\d+", colnames(A), value = TRUE)

# Select the corresponding columns in your data frame A
A_filtered <- A[, colnames(A) %in% selected_columns]

# Extract the "PI" values from the selected column names
pi_values <- paste0("PI",sub(".*_PI(\\d+)_.*", "\\1", selected_columns))

# Rename the columns in A_filtered with the cleaned "PI" values
colnames(A_filtered) <- pi_values

# Adding colnames
rownames(A_filtered) <- A[,1]
str(A_filtered)


# Set the threshold for proportion of zeroes you want to remove
threshold <- 0.5

# Calculate the proportion of zeroes in each column
zero_proportions <- colMeans(A_filtered == 0, na.rm = TRUE)

# Get the column indices to keep (where less than 50% are zeroes)
columns_to_keep <- which(zero_proportions < threshold)

# Subset A_filtered to include only the selected columns
A_filtered <- A_filtered[, columns_to_keep]





#write.csv(A_filtered,"result/A_filtered.csv", row.names = T)

##### Filtering based on lipids
lipid <- read.csv("data/lipid_class.csv")  

# Filtering only PC's
# Create a logical vector indicating rows where Compound_Name starts with "PC"
pc_names <- grepl("^PC", lipid$Compound_Name)

# Subset the data frame to include only rows with PC names
pc_data <- lipid[pc_names, ]

# View the first few rows of the subsetted data frame
str(pc_data)

unique_scans <- pc_data %>%
  group_by(X.Scan.) %>%
  summarize(Compound_Name = toString(unique(Compound_Name)))

# Create a function to add suffixes to duplicate row names
add_suffix <- function(names) {
  counts <- table(names)
  suffixes <- ave(seq_along(names), names, FUN = function(x) seq_len(length(x)) - 1)
  suffixes[suffixes > 0] <- sprintf("_%d", suffixes[suffixes > 0])
  paste(names, suffixes, sep = "")
}

# Apply the function to Compound_Name in unique_scans
unique_scans$Compound_Name <- add_suffix(unique_scans$Compound_Name)

# Subset rows in A_filtered using row names from unique_scans
subset_A <- A_filtered[rownames(A_filtered) %in% unique_scans$X.Scan., ]

# Replace row names with the updated Compound_Name
rownames(subset_A) <- unique_scans$Compound_Name


# Transpose the data frame A_PC
A_PC_transposed <- as.data.frame(t(subset_A))

# Identify rows with "CHECK" in row names and remove them
A_PC_transposed <- A_PC_transposed[!grepl("CHECK", rownames(A_PC_transposed)),, drop = FALSE]




# Filtering only PC's
# Create a logical vector indicating rows where Compound_Name starts with "PC"
lpc_names <- grepl("LPC", lipid$Compound_Name)

# Subset the data frame to include only rows with PC names
lpc_data <- lipid[lpc_names, ]

# View the first few rows of the subsetted data frame
str(lpc_data)

unique_scans <- lpc_data %>%
  group_by(X.Scan.) %>%
  summarize(Compound_Name = toString(unique(Compound_Name)))

# Create a function to add suffixes to duplicate row names
add_suffix <- function(names) {
  counts <- table(names)
  suffixes <- ave(seq_along(names), names, FUN = function(x) seq_len(length(x)) - 1)
  suffixes[suffixes > 0] <- sprintf("_%d", suffixes[suffixes > 0])
  paste(names, suffixes, sep = "")
}

# Apply the function to Compound_Name in unique_scans
unique_scans$Compound_Name <- add_suffix(unique_scans$Compound_Name)

# Subset rows in A_filtered using row names from unique_scans
subset_A <- A_filtered[rownames(A_filtered) %in% unique_scans$X.Scan., ]

# Replace row names with the updated Compound_Name
rownames(subset_A) <- unique_scans$Compound_Name


# Transpose the data frame A_lpc
A_LPC_transposed <- as.data.frame(t(subset_A))

# Identify rows with "CHECK" in row names and remove them
A_LPC_transposed <- A_LPC_transposed[!grepl("CHECK", rownames(A_LPC_transposed)),, drop = FALSE]


#Combining the data
A_PC_LPC <- cbind(A_PC_transposed,A_LPC_transposed)

write.csv(A_PC_LPC,"A_PC_LPC.csv", row.names = T)



