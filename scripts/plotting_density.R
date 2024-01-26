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
for (i in 6:10){
  d <- density(PC_A[,i])
  plot(d, main = colnames(PC_A[i]), col="red", ylim = c(0, max(d$y)))
  e <- density(log10_transformed_df[,i])
  plot(e, main=colnames(PC_A[i]), col="blue", ylim = c(0, max(e$y)))
  f <- density(centered_log10_df[,i])
  plot(f, main=colnames(PC_A[i]), col="green", ylim = c(0, max(f$y)))
}

par(mfrow = c(5, 3))

library(ggplot2)


library(sm)
attach(mtcars)

# create value labels
cyl.f <- factor(cyl, levels= c(4,6,8),
                labels = c("4 cylinder", "6 cylinder", "8 cylinder"))

# plot densities
sm.density.compare(mpg, cyl, xlab="Miles Per Gallon")
title(main="MPG Distribution by Car Cylinders")

# add legend via mouse click
colfill<-c(2:(2+length(levels(cyl.f))))
legend(locator(1), levels(cyl.f), fill=colfill)
