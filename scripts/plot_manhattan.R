#Package
library("CMplot")
library(vroom)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot)
library(data.table)

#Need LMM GWAS file all assoc
dir <- c("/Users/nirwantandukar/Documents/SAP genotype/gwas_results/LowP")

# List the files in the directory
file_list <- list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty list to store the data frames
df_list <- list()

# Read each file and store it in the list with the appropriate name
for (file in file_list) {
  # Extract the filename without the .txt extension
  file_name <- sub("\\.txt$", "", basename(file))
  
  # Read the file with vroom
  df <- vroom(file)
  
  # Subset columns 2, 1, 3, and 10
  #df <- df[, c(2, 1, 3, 10)]
  
  # Rename the columns
  colnames(df) <- c("SNP", "Chromosome", "Position", paste0(file_name))
  
  # Add the data frame to the list
  df_list[[file_name]] <- df
}


# Create PC1 data frame
PC1 <- df_list[[1]][, c(1:3, 4)]  # First sublist with columns 1-3 and 4
for (i in 2:4) { #for (i in 2:4) {
  PC1 <- cbind(PC1, df_list[[i]][, 4])  # Append the fourth column from the remaining sublists
  colnames(PC1)[ncol(PC1)] <- colnames(df_list[[i]])[4]  # Use the column name from the sublist
}
#PC1$SNP <- paste("SNP", PC1$Position, sep="_")
head(PC1)


# Create PC2 data frame
PC2 <- df_list[[5]][, c(1:3, 4)]  # Fifth sublist with columns 1-3 and 4
for (i in 6:8) {
  PC2 <- cbind(PC2, df_list[[i]][, 4])  # Append the fourth column from the remaining sublists
  colnames(PC2)[ncol(PC2)] <- colnames(df_list[[i]])[4]  # Use the column name from the sublist
}
#PC2$SNP <- paste("SNP", PC2$Position, sep="_")




# Create PC3 data frame
PC3 <- df_list[[9]][, c(1:3, 4)]  # Ninth sublist with columns 1-3 and 4
for (i in 10:11) {
  PC3 <- cbind(PC3, df_list[[i]][, 4])  # Append the fourth column from the remaining sublists
  colnames(PC3)[ncol(PC3)] <- colnames(df_list[[i]])[4]  # Use the column name from the sublist
}
#PC3$SNP <- paste("SNP", PC3$Position, sep="_")


# SNPs_1 <- PC1[
#   PC1$PC_15_0__0_0 < 1e-9 |
#     PC1$PC_16_0_20_5 < 1e-9 |
#     PC1$PC_16_1__16_1 < 1e-9 |
#     PC1$PC_16_1__20_4__PC_16_0__20_5__PC_18_2__18_3  < 1e-7 ,1
# ]



#quartz()
setwd("~/Desktop")

CMplot(PC1, plot.type="m",multracks=TRUE,threshold=c(1e-8,1e-6),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=0.5, file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)


CMplot(PC2, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=0.5, file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)

CMplot(PC3, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=0.5, file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)

