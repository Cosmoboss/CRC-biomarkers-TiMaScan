#Go to the directory which contains the files
setwd("C:/Users/Gebruiker/Documents/LUMC/secreted_mono_genes")
#Read in all the gene expression data from all cancers
TCGA <- read.table("TCGA_genes_after_filters.xls", sep="\t", header=T)
#Get all of the unique cancers that are listed in the TCGA data
cancers <- unique(TCGA$cancer)

#Loop through all the gene expression values per cancer
for (x in 1:length(cancers)){
  #Select all values of the current cancer
  c <- TCGA[which(TCGA$cancer == cancers[x]),- length(TCGA)]
  c <- t(c) #invert the table
  #Remove any rows that are NA
  row.has.na <- apply(c, 1, function(i){any(is.na(i))}) 
  c <-  c[!row.has.na,]
  #Calculate the means of each gene in the current cancer type
  means <- rowMeans(c)
  
  if (x == 1){ #Make a variable called all_means in the first time looping
    all_means <- means
  }else{ #Merge all the means of all genes together from all cancers
    all_means <- merge(all_means, means, by.x='row.names', by.y='row.names', all=TRUE)
    row.names(all_means) <- all_means[, 1]
    all_means <- all_means[,- 1]
  }
}
 #Change the colnames of the means to their corrisponding cancer type
colnames(all_means) <- cancers
#make a correlation heatmap of all the average gene expression values in all cancers
heatmap(cor(v, use = 'complete.obs'))
