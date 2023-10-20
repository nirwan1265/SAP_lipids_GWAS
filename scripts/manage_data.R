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
B_filtered <- B_filtered[, columns_to_keep_B]


#write.csv(A_filtered,"results/A_filtered.csv", row.names = T)
#write.csv(B_filtered,"results/B_filtered.csv", row.names = T)

##### Filtering based on lipids
lipid <- read.csv("data/lipid_class.csv")  
lipid$X.Scan. <- as.character(lipid$X.Scan.)
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

unique_scans_all_lipid <- lipid %>%
  dplyr::group_by(X.Scan.) %>%
  dplyr::summarize(Compound_Name = toString(unique(Compound_Name)))



# Process duplicates
unique_scans_all_lipid$Compound_Name <- with(unique_scans_all_lipid, ave(Compound_Name, Compound_Name, FUN = function(x) {
  if (length(x) > 1) {
    return(paste0(x, "_", 1:length(x)))
  } else {
    return(x)
  }
}))

# Convert X.Scan. column to character type
unique_scans_all_lipid$X.Scan. <- as.character(unique_scans_all_lipid$X.Scan.)

# Subset rows in A_filtered using row names from unique_scans
#subset_A <- A_filtered[rownames(A_filtered) %in% unique_scans$X.Scan., ]
# Convert rownames to a column
A_filtered <- A_filtered %>%
  mutate(rowname_col = rownames(A_filtered))

B_filtered <- B_filtered %>%
  mutate(rowname_col = rownames(B_filtered))

# Perform the inner join
subset_A_all_lipid <- inner_join(A_filtered, lipid, 
                                 by = c("rowname_col" = "X.Scan."))

subset_B_all_lipid <- inner_join(B_filtered, lipid, 
                                 by = c("rowname_col" = "X.Scan."))

subset_A_all_lipid <- subset_A_all_lipid[, c(ncol(subset_A_all_lipid), 1:(ncol(subset_A_all_lipid)-1))]

subset_B_all_lipid <- subset_B_all_lipid[, c(ncol(subset_B_all_lipid), 1:(ncol(subset_B_all_lipid)-1))]


write.csv(subset_A_all_lipid,"SetA_all_lipids.csv")
write.csv(subset_B_all_lipid,"SetB_all_lipids.csv")

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



