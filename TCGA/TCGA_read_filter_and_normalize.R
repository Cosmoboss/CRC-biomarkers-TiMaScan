#Import all necessary libraries
library(edgeR)

#List the short names of the cancer types
names <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", 
           "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD",
           "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
           "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

#Loop through the list of cancers
for (y in names){
  #Go to the directory for a specific cancer
  setwd(paste0("C:/Users/Gebruiker/Documents/LUMC/TCGA/TCGA_RSEM/", y))
  #Read in the raw table of data from the current cancer type
  raw_table <- read.table(paste0(y,".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"), sep="\t", header=T, stringsAsFactors = F)
  #This data has three rows for each sample in the dataset
  #The raw count, the scaled estimate and the transcript ID
  #We only need the raw counts
  
  #Safe the data in another data matrix
  data_matrix <- raw_table
  #We want the gene names as row names
  row.names(data_matrix) <- raw_table[,1]
  #The column names are the sample names
  #The first row is like a second layer of columns, but we don't need those
  #In addition, we only want the count data and nothing else.
  #Starting from position 2, we pick all columns with count data using a step size of 3
  data_matrix <- data_matrix[- 1,seq(2,length(data_matrix),3)]
  
  #Safe the row names
  rn <- data_matrix
  #Change all count values to numeric counts 
  data_matrix <- as.matrix(sapply(data_matrix, as.numeric))
  #Re-add the rownames
  row.names(data_matrix) <- row.names(rn)
  
  #Safe the column names
  columns <- colnames(data_matrix)
  #Make an empty list to put in the sample names that are from primary tumors
  pt_list <- c()
  
  #Make a counting number for later to add primary tumor column positions
  num=0
  #Loop through all samples from the current cancer
  for (x in columns){
    num = num+1 #Add a column position of 1
    if (substr(x, 14,15) == "01"){ #If the position is primary tumor i.e. 01 at position 14-15
      pt_list <- c(pt_list, num) #Add the position to the list of primary tumors
    }
  }
  
  #From the data_matrix only pick up all the primary tumor samples
  data_matrix <- data_matrix[, pt_list]
  myCPM <- cpm(data_matrix) #Calculate the CPM
  thresh <- myCPM > 2 #Any gene with a CPM below 2 is considered absent/low expressed
  n=ncol(data_matrix)/2 #50% of all samples in the cancer (number differs per cancer)
  #If the gene has a CPM above 2 in 50% or more of the cancers keep it.
  keep <- rowSums(thresh) >= n 
  summary(keep) #To show how many genes were kept.
  
  #Voom normalize the count data
  voomed <- voom(data_matrix[keep,])
  
  #Safe the filtered counts in a file
  write.table(voomed$E, paste0(y, "_datamatrix_filtered.xls"), quote=F, sep="\t")
}
