# Load the necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot) # Load cowplot for plot_grid

# Log10 transform each column separately
log10_transformed_df <- log10(PC_A)

# Calculate the median per column for log10-transformed data
median_log10 <- apply(log10_transformed_df, 2, median, na.rm = TRUE)

# Center each column by subtracting the median
centered_log10_df <- sapply(1:ncol(log10_transformed_df), function(col_index) {
  log10_transformed_df[, col_index] - median_log10[col_index]
})


par(mfrow = c(5, 3))

d <- density(PC_A[,5])
plot(d)
e <- density(log10_transformed_df[,5])
plot(e)
f <- density(centered_log10_df[,5])
plot(f)



