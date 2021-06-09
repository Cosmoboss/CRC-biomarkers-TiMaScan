#Safe all the cancer types in a list
cancers <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COADREAD", "DLBC", 
           "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD",
           "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "SARC", "SKCM",
           "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

#Safe all the bodily systems for each cancer (same order as above) in a list
systems <- c("Endocrine system", "Urinary system", "Breast invasive carcinoma", "Reproductive system", "Accessory organs GI tract", "Colorectal",
             "Immune system", "Gastrointestinal tract", "Central nervous system", "Gastrointestinal tract", "Urinary system", "Urinary system",
             "Urinary system", "Central nervous system", "Accessory organs GI tract", "Respiratory system", "Respiratory system", "Respiratory system",
             "Reproductive system", "Accessory organs GI tract", "Endocrine system", "Reproductive system", "Sarcoma", "Melanomas",
             "Gastrointestinal tract", "Reproductive system", "Endocrine system", "Immune system", "Reproductive system", "Reproductive system", "Melanomas")

#Loop through all cancers
for (i in 1:length(cancers)){
  if (i == 1){ #The first iteration
    print(cancers[i]) #Print the current cancer type which is running
    #Go to the directory where the files are stored
    setwd(paste0("C:/Users/Gebruiker/Documents/LUMC/TCGA/TCGA_RSEM/", cancers[i]))
    #Read in the filtered datamatrix
    TCGA_all <- read.table(paste0(cancers[i],"_datamatrix_filtered.xls"), sep="\t", header=T)
    #Make a pmeta with all TCGA samples, cancer types and bodily system
    supplement_all <- data.frame(sample=colnames(TCGA_all), cancer=cancers[i], system=systems[i])
  }else{ #Any other iteration than the first
    print(cancers[i]) #Print the current cancer type which is running
    #Go to the directory where the files are stored
    setwd(paste0("C:/Users/Gebruiker/Documents/LUMC/TCGA/TCGA_RSEM/", cancers[i]))
    #Read in the table with the filtered genes from the current cancer type
    table <- read.table(paste0(cancers[i],"_datamatrix_filtered.xls"), sep="\t", header=T)
    #Merge all the files by gene names into one big data frame
    TCGA_all <- merge(TCGA_all, table, by.x="row.names", by.y="row.names", all=T)
    row.names(TCGA_all) <- TCGA_all[, 1]
    TCGA_all <- TCGA_all[,- 1]
    #Merge all pmeta data into one big data frame
    s <- data.frame(sample=colnames(table), cancer=cancers[i], system=systems[i])
    supplement_all <- rbind(supplement_all, s)
  }
}


#Filter out the secreted proteins and macrophage/monocyte genes
setwd("C:/Users/Gebruiker/Documents/LUMC/secreted_mono_genes")

secretory <- read.table("Uniprot_Secretory_Exosomes_genes", sep="\t", header=T)
mmt1 <- read.table("Combined_sp_Markers_10X_gt_0.5", sep="\t", header=T)
mmt2 <- read.table("Mono_Macro_sp_Markers_SS2_gt_0.5", sep="\t", header=T)
MMT <- rbind(mmt1, mmt2)

#The row names in TCGA_all contain the gene name, an '|' and a number. 
#This seperates the gene names which then can be filtered for secr and macr genes
TCGA_gene_list <- c()
for (x in row.names(TCGA_all)){
  y <- strsplit(x, "|", fixed = T)[[1]]
  TCGA_gene_list <- c(TCGA_gene_list, y[1])
}

TCGA_filtered <- c()
TCGA_selected <- c()
#Checks for presence of secreted proteins, if the gene is a ?. Any genes that
#don't meat this criteria are selected
for (x in TCGA_gene_list){
  if (x %in% secretory$Gene.names){
    TCGA_filtered <- c(TCGA_filtered, x)
  }else if (x == "?"){
    print(x)
  }else{
    TCGA_selected <- c(TCGA_selected, x)
  }
}

#Select all secreted proteins and write the expression values of those in a seperate file
TCGA_secr_proteins <- TCGA_all[TCGA_gene_list %in% TCGA_filtered,]
write.table(TCGA_secr_proteins, file="TCGA_secreted_proteins.xls", sep="\t", quote=F)

`%notin%` <- Negate(`%in%`) # Create a negate of %in% operator

#Select all genes that aren't in macrophages/monocytes
TCGA_selected_genes <- TCGA_selected[TCGA_selected %notin% row.names(MMT)]
#Now select all expression values of the genes that weren't filtered out
TCGA_after_filters <- TCGA_all[which(TCGA_gene_list %in% TCGA_selected_genes),]

#Invert the rows and columns and add another column that lists cancer types for each sample
TCGA <- cbind(t(TCGA_after_filters), supplement_all$cancer)
#Change the column name to cancer instead of "X"
colnames(TCGA)[colnames(TCGA) == ""] <- "cancer"
write.table(TCGA, "TCGA_genes_after_filters.xls", quote=F, sep="\t")
#Write the pmeta in a seperate file
write.table(supplement_all, "TCGA_pmeta.xls", quote=F, sep="\t")
