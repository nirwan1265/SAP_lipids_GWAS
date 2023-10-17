# Libraries
library(dplyr)

# Loading data
A <- read.csv("data/SetA_lipid_FLO2020Control.csv")
B <- read.csv("data/SetB_lipid_FLO2022_lowP.csv")

###### Filtering just the Peak intensities 

# Find column names that match the pattern "S1_Run" at the end
selected_columns_A <- grep("ISTD", colnames(A), value = TRUE)
selected_columns_B <- grep("ISTD", colnames(B), value = TRUE)

# Select the corresponding columns in your data frame A
A_filtered <- A[, colnames(A) %in% selected_columns_A]
B_filtered <- B[, colnames(B) %in% selected_columns_B]

# Extract the "PI" values from the selected column names
pi_values_A <- paste0("PI",sub(".*_PI(\\d+)_.*", "\\1", selected_columns_A))
pi_values_B <- paste0("PI",sub(".*_PI(\\d+)_.*", "\\1", selected_columns_B))
names(A_filtered)
intersect(pi_values_A,pi_values_B)


##### Filtering based on lipids
lipid <- read.csv("data/lipid_class.csv")  
unique_scans_all_lipid <- lipid %>%
  group_by(X.Scan.) %>%
  summarize(Compoundx_Name = toString(unique(Compound_Name)))
add_suffix <- function(names) {
  counts <- table(names)
  suffixes <- ave(seq_along(names), names, FUN = function(x) seq_len(length(x)) - 1)
  suffixes[suffixes > 0] <- sprintf("_%d", suffixes[suffixes > 0])
  paste(names, suffixes, sep = "")
}

# Apply the function to Compound_Name in unique_scans
unique_scans_all_lipid$Compound_Name <- add_suffix(unique_scans_all_lipid$Compound_Name)
