
#Read in the spectral files and get out the genes in them
spectral <- read.table("TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv", sep="\t", header=T)
spectral_genes <- spectral$Gene[- 1722] #The last row is a total of all genes

#Read in all the proteingroup files that contain all the information for the proteomics datasets
PXD000089 <- read.table("pxd000089_proteinGroups.txt", sep="\t", header=F)
PXD001626 <- read.table("pxd001626_proteinGroups.txt", sep="\t", header=F)
PXD005693 <- read.table("pxd005693_proteinGroups.txt", sep="\t", fill=TRUE)

#Get the column containing all gene names and ignore the first row which are column names
PXD000089_genes <- PXD000089[- 1,7]
PXD001626_genes <- PXD001626[- 1,7]
PXD005693_genes <- PXD005693[- 1,7]

#Make an empty list to contain all the seperated gene names
list <- c()
for (x in PXD000089_genes){ #Loop trough all the gene names
  y <- strsplit(x, ";", fixed = T)[[1]] #Split each gene name for ";"
  list <- c(list, y[1]) #Merge the list with the obtained gene names
}
PXD000089_genes <- unique(list[!is.na(list)]) #Put all unique genes and not the NA's back into the protein dataset gene list

#Repeat the steps for both other protein datasets
list <- c()
for (x in PXD001626_genes){ #Loop trough all the gene names
  y <- strsplit(x, ";", fixed = T)[[1]] #Split each gene name for ";"
  list <- c(list, y[1]) #Merge the list with the obtained gene names
}
PXD001626_genes <- unique(list[!is.na(list)]) #Put all unique genes and not the NA's back into the protein dataset gene list


list <- c()
for (x in PXD005693_genes){ #Loop trough all the gene names
  y <- strsplit(x, ";", fixed = T)[[1]] #Split each gene name for ";"
  list <- c(list, y[1]) #Merge the list with the obtained gene names
}
PXD005693_genes <- unique(list[!is.na(list)]) #Put all unique genes and not the NA's back into the protein dataset gene list


