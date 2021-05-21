#import all the packages
library(edgeR)
library(DESeq2)
library(ggplot2)
library(biomaRt)

#Go to the directory where the .count files are stored
setwd("C:/Users/Gebruiker/Documents/LUMC/GSE_identifier_script/RNA-Seq/gse74369")

#Make a pmeta data matrix corrisponding to the sample files
#Do note that this code for the pmeta relies on the tumor and normal samples being
#clustered together in the directory based on tissue and not randomized.
#Two of the same loops to make the C/N list are commented. 
#The difference is that the "N" and "C" are switched.
#You can comment/uncomment the code depending on which position the samples are listed in the files
c <- 9
n <- 7
samples_list <- c()
for (x in 1:n){samples_list[x] <- paste0("N",x)}
for (x in 1:c){samples_list[length(samples_list)+1] <- paste0("C",x)}

#for (x in 1:c){samples_list[x] <- paste0("C",x)}
#for (x in 1:n){samples_list[length(samples_list)+1] <- paste0("N",x)}

ClassMetaData_list <- c()
for (x in 1:n){ClassMetaData_list[x] <- "normal"}
for (x in 1:c){ClassMetaData_list[length(ClassMetaData_list)+1] <- "cancer"}

#for (x in 1:c){ClassMetaData_list[x] <- "cancer"}
#for (x in 1:n){ClassMetaData_list[length(ClassMetaData_list)+1] <- "normal"}

#Look at all the .count files in the current directoy, save the names and remove the extention
countf = list.files(pattern = "\\.count$")
countf = gsub(pattern = "\\.count$", "", countf)
#Add the sample names to the pmeta
pmeta <- data.frame(samples=samples_list, ClassMetaData=ClassMetaData_list, ID=countf)

#Also add the filenames with the extentions
pmeta$countf = paste(pmeta$ID,".count", sep="")
#Make a Digital gene expression file containing all the count data for the current study
counts = readDGE(pmeta$countf)$counts
#A few extra rows got added that summarize the data
#Make a list to see which data are these extra rows for later removal
noint = rownames(counts) %in%
  c("__no_feature","__ambiguous","too_low_aQual",
    "not_aligned","__alignment_not_unique")
#Calculate the counts per million for the count data
cpms = cpm(counts)
#Filter out the absent/low expressed genes from the datasets
#Keep the genes where the CPM is above 1 in 3 or more of the samples and ignore the extra rows
keep = rowSums(cpms>1)>=3 & !noint
#Safe all the count values of the genes that are not absent/low expressed
counts = counts[keep,]

#Safe the counts and pmeta in an extra variable
counts_safe <- counts
pmeta_safe  <- pmeta

#Safe the pmeta before removing outliers
write.table(pmeta_safe, file = "gse74369_pmeta_before_outliers.xls", sep = "\t", quote=F)
#Remove any outlier samples from the counts and pmeta if needed
#Put in all the positions of the outlier samples in the counts and pmeta
counts <- counts_safe[,- c(7,8,9,18,19,20)]
pmeta <- pmeta_safe[- c(7,8,9,18,19,20),]

#Safe the counts in a matrix
write.table(counts, file = "gse74369_counts.xls", sep = "\t", quote=F)

#Change the column names to the same sample names as in pmeta
colnames(counts) <- pmeta$samples

#Make a DGElist with the filtered count data
dge <- DGEList(counts = counts, group = pmeta$ClassMetaData)

#calculate the normalization factors
norm <- calcNormFactors(dge)

#Normalize the count data using voom
d_voom <- voom(norm, plot = TRUE)

#Make a PCA plot of the count data to look for potential outliers
exprset22=t(d_voom$E)
PC<-prcomp(exprset22)
#PCA plot
PCi<-data.frame(PC$x,Species=pmeta$ClassMetaData)
ggplot(PCi,aes(x=PC1,y=PC2,col=pmeta$ClassMetaData))+
  geom_point(size=3)+ #Dot size
  geom_text(size=5, aes(label=pmeta$samples),hjust = 0, nudge_x = 0.1) + #Text size
  stat_ellipse() + #Confidence interval ellipse at 100%
  stat_ellipse(level=0.80,linetype = 2)+ #Condifence interval ellipse at 80%
  scale_color_manual(values = c("blue", "black"), labels = c("Cancer", "Normal"))+ #Colours and cohorts
  labs(col = "Disease Status") + #Legend title
  theme(panel.background = element_blank(), #Extra bolding and positioning of x/y-axis
        axis.line  = element_line(colour = "black"),
        legend.text = element_text(face = "bold"),
        legend.key = element_rect(fill = 'white', size = 0.5, linetype='solid'),
        legend.key.size = unit(1.5, 'lines'),
        #legend.title=element_blank(),
        axis.text.x=element_text(size=rel(1.5), face="bold", colour = "black", angle = 30, hjust=1),
        axis.text.y=element_text(size=rel(1.5), face="bold", colour = "black"),
        axis.title=element_text(size=rel(1.5)))

#If there are no more outliers left
#Add gene names to the ensembl ID's
ensembl = useMart('ensembl', dataset = "hsapiens_gene_ensembl")
ensembl_gene_names <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), filters = "ensembl_gene_id", values = row.names(counts), mart= ensembl)
#Add the ensembl gene names to the data. Fill any row as NA if there was no gene names
rna_seq_processed <- merge(d_voom$E, ensembl_gene_names, by.x = 'row.names', by.y = "ensembl_gene_id", all = TRUE)
row.names(rna_seq_processed) <- rna_seq_processed[,1]
rna_seq_processed <- rna_seq_processed[- c(1:3),- 1]
write.table(rna_seq_processed, file = "gse74369_processed.xls", sep = "\t", quote=F)
write.table(pmeta, file = "gse74369_pmeta_after_outliers.xls", sep = "\t", quote=F)
