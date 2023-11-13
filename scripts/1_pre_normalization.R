library(stringr)
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


# Final Table
A <- cbind(A_id,A_peaks)

#Removing columns with InjBL cause SERRF does not want it
A <- A %>%
  dplyr::select(-contains("InjBL"))


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
num_cols <- ncol(A_peaks) - 1

# Create a new row with increasing numbers from 1 to num_cols
new_row <- c("time", 1:num_cols)

# Add the new row as the first row in A_peaks
A <- rbind(new_row, A)


write.csv(A,"SERRF_control.csv")
