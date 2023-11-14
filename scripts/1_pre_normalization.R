library(stringr)
library(purrr)
################################################################################
############ Loading the Quantitative Peaks Intensities.  ######################
################################################################################

# Control
A <- read.csv("data/SetA_lipid_FLO2019Control.csv")
# Low P
#B <- read.csv("data/SetB_lipid_FLO2022_lowP.csv")
# Separating the tables
A_peaks <- A[14:882]
A_id <- A[1]

# Extract everything including "S1," "InjBL," and "QC"
extracted_strings <- str_extract(names(A_peaks), "(S1_[^.]|InjBL.[^.]|QC_R[^.]|ISTD_[^.]+).*")

# Remove everything after ".mzML," including ".mzML.Peak.area"
extracted_strings <- sub("\\.mzML.*", "", extracted_strings)

# Print the extracted strings
head(extracted_strings)

# Replace column names
colnames(A_peaks) <- extracted_strings

# Extract "Run" and numbers from the updated column names
run_numbers <- as.numeric(sub(".*Run(\\d+).*", "\\1", extracted_strings))
head(run_numbers)

# Sort the extracted numbers and keep track of the original indices
sorted_indices <- order(run_numbers)
head(sorted_indices)
str(sorted_indices)

# Arrange the columns in A_peaks based on the sorted indices
A_peaks <- A_peaks[, sorted_indices]
A_peaks[1:5,1:10]

head(A_peaks)

# Final Table
A <- cbind(A_id,A_peaks)

# Removing compounds that have values 10 times less than blanks
# Calculate the minimum value for the "S1_Run*" columns only
A$Min_S1_Run <- apply(A[, grepl("^S1_Run", colnames(A))], 1, min)
# Calculate the average of InjBL.334_Run* columns
A$InjBL.334_Avg <- apply(A[, grepl("^InjBL\\.334_Run", colnames(A))], 1, max)
#A$InjBL.334_Avg <- rowMeans(A[, grepl("^InjBL\\.334_Run", colnames(A))])
# Create a logical vector to identify rows that meet the condition
condition_met <- A$Min_S1_Run >= 10 * A$InjBL.334_Avg

# Subset the data frame to keep only the rows where the condition is met
A <- A[condition_met, ]
# Remove the "Min_S1_Run" and "InjBL.334_Avg" columns
A_peaks <- A_peaks[, !grepl("^(Min_S1_Run|InjBL\\.334_Avg)$", colnames(A_peaks))]


#Removing columns with InjBL cause SERRF does not want it
A <- A %>%
  dplyr::select(-contains("InjBL"))


# Filtering relevant columns
relevant_cols <- grep("^S1_Run", names(A), value = TRUE)
A <- A %>%
  dplyr::select(all_of(relevant_cols)) %>%
  dplyr::select(where(~ mean(.x == 0) <= 0.3))


# Create a new row with labels based on column names
new_row <- data.frame("A")


# Add QC, Sample, or Validate labels based on column names
col_names <- colnames(A)
for (col_name in col_names) {
  if (grepl("QC_", col_name)) {
    new_label <- "qc"
  } else if (grepl("S1_", col_name)) {
    new_label <- "sample"
  } else if (grepl("ISTD_", col_name)) {
    new_label <- "validate"
  } else if (grepl("row.ID", col_name)) {
    new_label <- "label"
  } else {
    new_label <- ""
  }
  
  new_row[[col_name]] <- new_label
  
}
new_row <- as.data.frame(new_row[,-1])

# Add the new row to the data frame
A <- rbind(new_row, A)

# Adding sample numbers
# Get the number of columns in A_peaks
num_cols <- ncol(A) - 1

# Create a new row with increasing numbers from 1 to num_cols
new_row <- c("time", 1:num_cols)

# Add the new row as the first row in A_peaks
A <- rbind(new_row, A)


