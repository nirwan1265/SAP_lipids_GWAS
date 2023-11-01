# Load necessary libraries
library(dplyr)
library(vroom)

# Define the directory path and file pattern
directory_path <- "/Users/nirwantandukar/Documents/SAP genotype/gwas_results/LowP"
file_pattern <- ".txt"

# List the files matching the pattern
file_list <- list.files(path = directory_path, pattern = file_pattern)

# Loop through the files, read and transform them, then assign to variables
for (file_name in file_list) {
  # Remove the ".txt" extension from the variable name
  variable_name <- sub("\\.txt$", "", file_name)
  
  # Read the file content using vroom
  data <- vroom(file.path(directory_path, file_name))
  
  # Transform the data
  transformed_data <- data %>%
    mutate(rs = paste0("SNP_", ps)) %>%
    select(Marker = rs, Chr = chr, Pos = ps, P = p_lrt)
  
  # Assign the transformed data to a variable with the cleaned name
  assign(variable_name, transformed_data)
}

# Display transformed data for one of the files
head(PC_16_1_18_1)



getwd()
# Saving
# doing individually cause R crashes
write.table(PC_16_0_18_1,"PC_16_0_18_1_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_16_0_18_1_2_3_hexadecanoyloxy_2_octadec_9_enoyloxy_propyl_phosphonato_oxy_ethyl_trimethylazanium,"PC_16_0_18_1_2_3_hexadecanoyloxy_2_octadec_9_enoyloxy_propyl_phosphonato_oxy_ethyl_trimethylazanium_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_16_1_18_1,"PC_16_1_18_1_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_0_18_1,"PC_18_0_18_1_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_1_18_1,"PC_18_1_18_1_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_1_18_2,"PC_18_1_18_2_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_1_18_2_trimethyl_2_3_octadec_11_enoyloxy_2_octadeca_9_12_dienoyloxy_propyl_phosphonato_oxy_ethyl_azanium,"PC_18_1_18_2_trimethyl_2_3_octadec_11_enoyloxy_2_octadeca_9_12_dienoyloxy_propyl_phosphonato_oxy_ethyl_azanium_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_1_20_1,"PC_18_1_20_1_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_1_20_3,"PC_18_1_20_3_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_1_20_4,"PC_18_1_20_4_filtered.txt",quote=F, row.names=F, sep=" ")
write.table(PC_18_1_24_1,"PC_18_1_24_1_filtered.txt",quote=F, row.names=F, sep=" ")

