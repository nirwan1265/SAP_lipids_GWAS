#!/bin/bash

# Sample data file
data_file="your_data_file.txt"

# Loop through chromosomes from 1 to 10
for chr in {1..10}; do
  # Define the output file name with leading zeros
  output_file="PC_18_2__0_$(printf "%03d" $chr).txt"

  # Extract and transform the data for the current chromosome
  awk -v chr="$chr" 'BEGIN {print "Marker Chr Pos p"} $1 == chr { print "SNP_"$3, $1, $3, $10 }' "$data_file" > "$output_file"

  echo "Processed chromosome $chr and saved to $output_file"
done