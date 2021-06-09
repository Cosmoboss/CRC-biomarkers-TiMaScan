#Import all necessary libraries
library(limma)

#Go to the directory which has the file with all filtered TCGA genes/cancers
setwd("C:/Users/Gebruiker/Documents/LUMC/secreted_mono_genes")
TCGA <- read.table("TCGA_genes_after_filters.xls", sep="\t", header=T)

#Get all unique cancer types in a lsit
cancers <- unique(TCGA$cancer)

#Go to the directory where all the files will be safed
setwd("C:/Users/Gebruiker/Documents/LUMC/TCGA/TCGA_vs_COADREAD_topTables")

#Get the expression values from genes in  from colorectal cancer
COADREAD <- TCGA[which(TCGA$cancer == "COADREAD"),- length(TCGA)]
COADREAD <- t(COADREAD) #invert the table
#Check for any gene that isn't present in CRC
row.has.na <- apply(COADREAD, 1, function(x){any(is.na(x))})
#Remove genes that aren't present in CRC
COADREAD <-  COADREAD[!row.has.na,]
#Remove COADREAD from the cancers list
cancers <- cancers[- which(cancers == "COADREAD")]

#loop through all cancers other than CRC
for (x in 1:length(cancers)){
  #fetch the expression values from the current cancer
  c <- TCGA[which(TCGA$cancer == cancers[x]),- length(TCGA)]
  c <- t(c)#invert the table
  
  #Check if there's any outliers present and remove them
  row.has.na <- apply(c, 1, function(x){any(is.na(x))})
  c <-  c[!row.has.na,]
  
  #Fetch all the genes that are present in both types of cancers
  genes <- c[which(row.names(c) %in% row.names(COADREAD)),]
  #Merge the two cancer types into one big data frame
  merged <- merge(COADREAD, genes, by.x='row.names', by.y='row.names')
  row.names(merged) <- merged[, 1]
  merged <- merged[,- 1]
  
  #Make a model
  CancerType_list <- c()
  for (i in 1:dim(COADREAD)[2]){CancerType_list[i] <- "COADREAD"}
  for (i in 1:dim(c)[2]){CancerType_list[length(CancerType_list)+1] <- "c"}
  samples=as.factor(CancerType_list)
  mod<-model.matrix(~0+samples)
  
  #Do differential gene expression analysis between CRC and the other cancer
  fit <- lmFit(merged, mod)
  #change the sample names according to the model.matrix
  contr <- makeContrasts(samplesCOADREAD - samplesc, levels = colnames(coef(fit)))
  fit <- contrasts.fit(fit, contr)
  fit <- eBayes(fit)
  tT <- topTable(fit,number = Inf, p.value=0.05, lfc=0.5, adjust.method = "fdr")
  
  #Calculate the mean expression value for each cancer type
  COADREAD_means <- rowMeans(COADREAD)
  c_means <- rowMeans(genes)
  
  #Calculate SD in both cancer types
  COADREAD_SD=apply(as.matrix(COADREAD),1, sd)
  c_SD=apply(as.matrix(genes),1, sd)
  
  #Calculate variance in both cancer types
  COADREAD_vars <- rowSums((COADREAD - rowMeans(COADREAD))^2)/(dim(COADREAD)[2] - 1)
  c_vars <- rowSums((genes - rowMeans(genes))^2)/(dim(genes)[2] - 1)
  
  #Combine the toptable with the means, SD and variance of both cancers
  crc <- cbind(COADREAD_means, COADREAD_SD, COADREAD_vars)
  other <- cbind(c_means, c_SD, c_vars)
  crc_other <- merge(crc, other, by.x='row.names', by.y='row.names')
  row.names(crc_other) <- crc_other$Row.names
  crc_other <- crc_other[,- 1]
  tT_crc_other <- merge(tT, crc_other, by.x = 'row.names', by.y ='row.names')
  row.names(tT_crc_other) <- tT_crc_other$Row.names
  tT_crc_other <- tT_crc_other[,- 1]
  
  #Change the column names to something more logical and fitting to the cancers
  names <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "COADREAD_means", "COADREAD_SD", "COADREAD_vars",
             paste0(cancers[x], "_means"), paste0(cancers[x], "_SD"), paste0(cancers[x], "_vars"))
  colnames(tT_crc_other) <- names
  
  #Safe the toptable in a file
  write.table(tT_crc_other, paste0(cancers[x], "_vs_COADREAD_topTable.xls"), quote=F, sep = "\t")
}

#To check which genes are upregulated in CRC compared to other types of cancer
for (x in 1:length(cancers)){
  #Read in the toptable for each cancer
  toptable <- read.table(paste0(cancers[x], "_vs_COADREAD_topTable.xls"), header=T, sep = "\t")
  toptable <- toptable[toptable$P.Value < 0.05,] #select for pvalue below 0.05
  toptable <- toptable[toptable$logFC > 0.5,] #Select for logFC above 0.5
  #Add cancer types as colnames for easier merging of columns later on
  colnames(toptable) <- c(cancers[x], "b", "c","d")
  if (x == 1){
    upregulated_all <- toptable[1] #Add all upregulated genes to a dataframe
  }else{
    upregulated_all <- merge(upregulated_all, toptable[1], by="row.names", all=TRUE)
    row.names(upregulated_all) <- upregulated_all[,1]
    upregulated_all <- upregulated_all[,- 1]
  }
}

#Make a presence absence matrix with zero's and one's
u_all <- upregulated_all
u_all[!is.na(u_all)] <- 1
u_all[is.na(u_all)] <- 0

write.table(u_all, "all_cancers_1_vs_0_matrix.xls", quote=F, sep="\t")

