library(stringr)
library(purrr)
################################################################################
############ Loading the Quantitative Peaks Intensities.  ######################
################################################################################

# Low P
B <- read.csv("data/SetB_lipid_FLO2022_lowP.csv")
# Separating the tables
B_peaks <- B[14:791]
B_id <- B[1]

# Extract everything including "S1," "InjBL," and "QC"
extracted_strings <- str_extract(names(B_peaks), "(S1_[^.]|InjBL.[^.]|QC_R[^.]|ISTD_[^.]+).*")

# Remove everything after ".mzML," including ".mzML.Peak.area"
extracted_strings <- sub("\\.mzML.*", "", extracted_strings)

# Print the extracted strings
head(extracted_strings)

# Replace column names
colnames(B_peaks) <- extracted_strings

# Extract "Run" and numbers from the updated column names
run_numbers <- as.numeric(sub(".*Run(\\d+).*", "\\1", extracted_strings))
head(run_numbers)

# Sort the extracted numbers and keep track of the original indices
sorted_indices <- order(run_numbers)
head(sorted_indices)
str(sorted_indices)

# Arrange the columns in B_peaks based on the sorted indices
B_peaks <- B_peaks[, sorted_indices]
B_peaks[1:5,1:10]

head(B_peaks)

# Final Table
B <- cbind(B_id,B_peaks)
names(B)

# Removing compounds that have values 10 times less than blanks
# Calculate the minimum value for the "S1_Run*" columns only
B$Min_S1_Run <- apply(B[, grepl("^S1_Run", colnames(B))], 1, min)
# Calculate the average of InjBL.334_Run* columns
B$InjBL.334_Avg <- apply(B[, grepl("^InjBL\\.334_Run", colnames(B))], 1, max)
#B$InjBL.334_Avg <- rowMeans(B[, grepl("^InjBL\\.334_Run", colnames(B))])

# Create B logical vector to identify rows that meet the condition
condition_met <- B$Min_S1_Run >= 10 * B$InjBL.334_Avg

# Subset the data frame to keep only the rows where the condition is met
B <- B[condition_met, ]

# Remove the "Min_S1_Run" and "InjBL.334_Avg" columns
B <- B[, !grepl("^(Min_S1_Run|InjBL\\.334_Avg)$", colnames(B))]


#Removing columns with InjBL cause SERRF does not want it
B <- B %>%
  dplyr::select(-contains("InjBL"))


# Filtering relevant columns
# relevant_cols <- grep("^S1_Run", names(B), value = TRUE)
# B <- B %>%
#   dplyr::select(all_of(relevant_cols)) %>%
#   dplyr::select(where(~ mean(.x == 0) <= 0.3))


# Create B new row with labels based on column names
new_row <- data.frame("B")


# Add QC, Sample, or Validate labels based on column names
col_names <- colnames(B)
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
B <- rbind(new_row, B)

# Adding sample numbers
# Get the number of columns in B_peaks
num_cols <- ncol(B) - 1

# Create B new row with increasing numbers from 1 to num_cols
new_row <- c("time", 1:num_cols)

# Add the new row as the first row in B_peaks
B <- rbind(new_row, B)

# Saving
write.csv(B,"results/SERRF/Control_SEERF.csv")

# Need some post processing for the SERRF
#batch	B	B	B	B	B	B	B	B	B	B
#sampleType	qc	validate	sample	sample	sample	sample	sample	sample	sample	sample
#time	1	2	3	4	5	6	7	8	9	10
#No	label	QC000	sample01	GB001617	GB001333	GB001191	GB001827	GB001722	GB001468	GB001543	GB001347
#1	1_ISTD Ceramide (d18:1/17:0) [M+HCOO]- 	167879	185671	158256	164492	155000	150957	134195	184272	165878	157758

# Website for running SERRF:
#https://slfan.shinyapps.io/ShinySERRF/
