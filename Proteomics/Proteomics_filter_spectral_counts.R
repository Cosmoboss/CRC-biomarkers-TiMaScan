#import the libraries that will be used
library(edgeR)

#Go to the directory with the spectral coutn files
setwd("C:/Users/Gebruiker/Documents/LUMC/Spectral_file")

#Read in the spectral count files
spectral <- read.table("TCGA_Colon_VU_Proteome_CDAP.r2.spectral_counts.tsv", header = T, sep="\t")

#Safe the spectral data matrix in another variable
table <- spectral
#Add the gene names as row.names to the new variable
row.names(table) <- spectral[,1]
#Remove the last row (which is a total row) and the last 6 rows which don't contain count data
table <- table[- 1722,- c(1, 192:197)]
#Pick up all the unshared columns.
unshared <- table[,seq(2,190,2)]
#Calculate the CPMs for the unshared count table
cpms <- cpm(unshared)
#Keep any genes of which the CPMs are above 0,05 in at least 25% of the samples
keep <- rowSums(cpms>0.05)>23.75
filtered <- unshared[keep,]
#Make a voom plot with the unfiltered counts
voom(unshared, plot = T)
#Make a voom plot with the filtered counts
voom(filtered, plot = TRUE)

