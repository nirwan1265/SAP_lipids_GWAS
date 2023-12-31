# Subsetting just the genename from the gff3 files

split.names <- function(x,split){
  if(organism == "Sorghum"){
    split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
    split.genename2 <- unlist(strsplit(split.genename, split = ":", fixed = TRUE))[2]
  }
  else {
    split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
    split.genename2 <- unlist(strsplit(split.genename, split = "=", fixed = TRUE))[2]
  }
  return(split.genename2)
}


# Converting pvalues to Z-scores
zval <- function(x, output){
  pvalue <- unlist(as.numeric(x[4]))
  o <- p.to.Z(pvalue)
  return(o)
}
preprocess <- function(path, phenoname, n, organism){
  #Helper lists 
  zstat_list <- list()
  pvalue_list <- list()
  marker_list <- list()
  
  # Set working directory
  setwd(path)
  #file_path = paste0(path,"/GWAS/sorghum/")
  #file_path = paste0(path,"/GWAS/maize/")
  # Preprocess for Sorghum
  if(organism == "Sorghum"){
    file_path = paste0(path,"/GWAS/sorghum/")
    # Load Sorghum Genomic Ranges
    a <- 1
    for(i in sprintf("%02d",1:chr)){
      assign(paste0("gr.db",a), readRDS(paste0("GenomicRanges/sorghum/gr.db",i,".RDS")))
      a = a + 1
    }
    
    # Data processing before GBJ
    a <- 1
    file_list <- list.files(path = file_path, pattern = phenoname)
    for(i in 1:chr){#length(file_list)){
      assign(file_list[i], vroom(paste0(file_path,file_list[i])))
      d = as.data.frame(get(file_list[i]))
      d$zstat = unlist(apply(d,1,zval))
      names(d) <- c("Marker","chr","Start_Position","pvalue","Zvalue")
      d <- d %>% mutate_at(c('chr','Start_Position','pvalue','Zvalue'),as.numeric)
      d <- as.data.frame(d)
      assign(paste0("gr.q", i) , GenomicRanges::GRanges(seqnames = paste0("chr",sprintf("%02d",i)), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
      a = a + 1
      assign(paste0("common",i), as.data.frame(IRanges::findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
      
      e = get(paste0("common",i))
      e = e[which(e$first.X.Region == "gene"), ]
      e = e[,c(7,16,17,18)]
      colnames(e) = c("Gene","Marker","pvalue","zstat")
      assign(paste0("filter_common", i), e)
      assign(paste0("zstat",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
      assign(paste0("Marker",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
      assign(paste0("pvalue",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
      assign(paste0("genename",i),apply(get(paste0("Marker",i)),1,split.names))
      f <- get(paste0("zstat",i))
      f[,1]<- get(paste0("genename",i))
      f <- as.data.frame(t(f))
      colnames(f) <- f[1,]
      f <- f[-1,]
      f <- f %>% mutate_if(is.character,as.numeric, na.rm = T)
      f <- f[mixedsort(row.names(f)), ]
      assign(paste0("zstat",i),f)
      g <- get(paste0("Marker",i))
      g[,1]<- get(paste0("genename",i))
      g <- as.data.frame(t(g))
      colnames(g) <- g[1,]
      g <- g[-1,]
      g <- g[mixedsort(row.names(g)), ]
      assign(paste0("Marker",i),g)
      #return (get(paste0("Marker",i)))
      h <- get(paste0("pvalue",i))
      h[,1]<- get(paste0("genename",i))
      h <- as.data.frame(t(h))
      colnames(h) <- h[1,]
      h <- h[-1,]
      h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
      h <- h[mixedsort(row.names(h)), ]
      assign(paste0("pvalue",i),h)
    }
  }
  else if(organism == "Zea"){
    file_path = paste0(path,"/GWAS/maize/")
    # Load maize Genomic Ranges
    a <- 1
    for(i in sprintf("%02d",1:chr)){
      assign(paste0("gr.db",a), readRDS(paste0("GenomicRanges/maize/gr.db",i,".RDS")))
      a = a + 1
    }
    
    # Data processing before GBJ
    a <- 1
    file_list <- list.files(path = file_path, pattern = phenoname)
    for(i in 1:chr){#length(file_list)){
      assign(file_list[i], vroom(paste0(file_path,file_list[i])))
      d = as.data.frame(get(file_list[i]))
      d$zstat = unlist(apply(d,1,zval))
      names(d) <- c("Marker","chr","Start_Position","pvalue","Zvalue")
      d <- d %>% mutate_at(c('chr','Start_Position','pvalue','Zvalue'),as.numeric)
      d <- as.data.frame(d)
      assign(paste0("gr.q", i) ,GenomicRanges::GRanges(seqnames = paste0("chr",sprintf("%02d",i)), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
      a = a + 1
      assign(paste0("common",i), as.data.frame(IRanges::findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
      e = get(paste0("common",i))
      e = e[which(e$first.X.Region == "gene"), ]
      e = e[,c(7,16,17,18)]
      colnames(e) = c("Gene","Marker","pvalue","zstat")
      assign(paste0("filter_common", i), e)
      assign(paste0("zstat",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
      assign(paste0("Marker",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
      assign(paste0("pvalue",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
      assign(paste0("genename",i),apply(get(paste0("Marker",i)),1,split.names))
      f <- get(paste0("zstat",i))
      f[,1]<- get(paste0("genename",i))
      f <- as.data.frame(t(f))
      colnames(f) <- f[1,]
      f <- f[-1,]
      f <- f %>% mutate_if(is.character,as.numeric, na.rm = T)
      f <- f[mixedsort(row.names(f)), ]
      assign(paste0("zstat",i),f)
      g <- get(paste0("Marker",i))
      g[,1]<- get(paste0("genename",i))
      g <- as.data.frame(t(g))
      colnames(g) <- g[1,]
      g <- g[-1,]
      g <- g[mixedsort(row.names(g)), ]
      assign(paste0("Marker",i),g)
      #return (get(paste0("Marker",i)))
      h <- get(paste0("pvalue",i))
      h[,1]<- get(paste0("genename",i))
      h <- as.data.frame(t(h))
      colnames(h) <- h[1,]
      h <- h[-1,]
      h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
      h <- h[mixedsort(row.names(h)), ]
      assign(paste0("pvalue",i),h)
    }
  }
  zstat_list = list()
  for (i in 1:chr){
    zstat_list <- c(zstat_list, list(get(paste0("zstat",i))))
    names(zstat_list)[i] <- paste0("zstat",i)
  }
  for (i in 1:chr){
    pvalue_list <- c(pvalue_list, list(get(paste0("pvalue",i))))
    names(pvalue_list)[i] <- paste0("pvalue",i)
  }
  for (i in 1:chr){
    marker_list <- c(marker_list, list(get(paste0("Marker",i))))
    names(marker_list)[i] <- paste0("markers",i)
  }
  return (list(Zstat = zstat_list, pvalue = pvalue_list,Marker = marker_list))
}

path <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
phenoname="PC_18_1_20_1_"
organism <- "Sorghum"
chr <- 10
preprocess_data <- preprocess(path, phenoname, chr,  organism)

head(str(preprocess_data$pvalue))



pvalue_list <- preprocess_data[["pvalue"]]


# Create an empty data frame to store the results
result_df <- data.frame(Gene = character(0), pvalue = numeric(0))

# Loop through the list of p-value data frames
for (i in 1:length(pvalue_list)) {
  pvalue_df <- pvalue_list[[i]]
  
  # Calculate the minimum p-value for each gene
  min_pvalues <- apply(pvalue_df, 2, function(col) min(col, na.rm = TRUE))
  
  # Combine gene names and their corresponding minimum p-values into a data frame
  gene_pvalue_df <- data.frame(Gene = names(min_pvalues), pvalue = min_pvalues)
  
  # Sort the data frame by minimum p-values in ascending order
  gene_pvalue_df <- gene_pvalue_df[order(gene_pvalue_df$pvalue), ]
  
  # Get the top 2 genes with the lowest p-values
  top_genes <- head(gene_pvalue_df, 2)
  
  # Add the results to the final data frame
  result_df <- rbind(result_df, top_genes)
}

# Print the final result
print(result_df)



# Assuming your data frame is called preprocess_data
# Transpose the data frame
transposed_data <- as.data.frame(t(preprocess_data[["pvalue"]][["pvalue7"]]))
str(transposed_data)
head(transposed_data)


# Function to get the column name with the minimum value in each row
get_min_col <- function(row) {
  return(names(row)[which.min(row)])
}

# Apply the function row-wise to find the column with the minimum value in each row
min_cols <- apply(transposed_data, 1, get_min_col)

# Reorder the columns based on the minimum column for each row
transposed_data <- transposed_data[, min_cols]

# Print the first few rows to verify the result
head(transposed_data)

str(transposed_data)
transposed_data[1:5,1:5]


# Create an empty data frame to store the result
result_df <- data.frame(Gene = rownames(transposed_data), pvalue = NA)

# Loop through each row
for (i in 1:nrow(transposed_data)) {
  # Find the column name with the minimum value in the current row
  min_col <- names(which.min(transposed_data[i, ]))
  
  # Assign the minimum pvalue to the result data frame
  result_df$pvalue[i] <- transposed_data[i, min_col]
}

# Print the result data frame
print(result_df)
