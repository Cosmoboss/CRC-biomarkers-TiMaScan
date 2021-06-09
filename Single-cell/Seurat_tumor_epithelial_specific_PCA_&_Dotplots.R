#Import the necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)


#Load in the gse_integrated seurat object containing the anchor merged studies
load("gse_integrated.Rda")


# This code has several plots that are made using the integrated seurat object of the three studies
# In order, the code 
#
# 1. PCA of the average expression of each study cell type
#
# 2. Dotplot of gene expression in epithelial cells
#
# 3. Dotplot of tumor epithelial specific genes in all cells

#The study cell types need to be reformated first for a better visualization of the dotplots

my_levels = c("FRZ Tumor_Epithelial cells", "KUL Tumor_Epithelial cells", "SMC Tumor_Epithelial cells",
              "FRZ Normal_Epithelial cells", "KUL Normal_Epithelial cells", "SMC Normal_Epithelial cells",
              "FRZ Tumor_Stromal cells", "KUL Tumor_Stromal cells", "SMC Tumor_Stromal cells",
              "FRZ Normal_Stromal cells", "KUL Normal_Stromal cells", "SMC Normal_Stromal cells")
Idents(gse_integrated) <- factor(Idents(gse_integrated), levels= my_levels)

#===================== PCA of the average expression of each study cell type ==================#

#obtain the cluster averages expression of each identity class for all genes
cluster_averages <- AverageExpression(gse_integrated, group.by = "study_cell_type")

#Both the raw merged and integrated merged counts are given, therefore we select the second variable
#Which is the integrated table
integrated <- cluster_averages[[2]]
#Construct a PCA plot of the study cell types
exprset22=t(integrated)
species=factor(colnames(integrated))
PC<-prcomp(exprset22)
PCi<-data.frame(PC$x,Species=colnames(integrated))
#Shapes and colours for the PCA
shapes <- c(15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17)
kleurtjes <- c("medium orchid", "blue", "dark orange", "forest green",
               "medium orchid", "blue", "dark orange", "forest green",
               "medium orchid", "blue", "dark orange", "forest green")


ggplot(PCi,aes(x=PC1,y=PC2,col=species))+
  geom_point(size=3,alpha=2, shape=shapes)+ #Size and alpha just for fun
  scale_color_manual(values = kleurtjes)+#, labels = c("Cancer", "Normal"))+ #your colors here
  labs(col = "Cohort") +
  #scale_shape_manual(values=c(3, 16, 17))+
  theme(panel.background = element_blank(),
        axis.line  = element_line(colour = "black"),
        legend.text = element_text(face = "bold"),
        legend.key = element_rect(fill = 'white', size = 0.5, linetype='solid'),
        legend.key.size = unit(1.5, 'lines'),
        #legend.title=element_blank(),
        axis.text.x=element_text(size=rel(1.5), face="bold", colour = "black", angle = 30, hjust=1),
        axis.text.y=element_text(size=rel(1.5), face="bold", colour = "black"),
        axis.title=element_text(size=rel(1.5)))


#================================== Dotplot of gene expression in epithelial cells ===============================================#

#First we will need all the CRC specific and upregulated genes in a list
#Safe all CRC specific genes in a list
crc_specific_genes <- c("APOBEC1", "ATOH1", "BTNL3", "CASP5", "CEL", "SCARNA17",
                        "CHST5", "NOXO1", "PHGR1","TM4SF20", "ZNF573", "SLC6A4",
                        "TBX10", "MOGAT2", "NOX1", "SLC26A3", "TUBAL3", "EFNA2",
                        "GUCY2C", "MEP1A", "LY6G6D", "TLX1", "PMFBP1", "FBXO48",
                        "ZNF699", "DUSP27", "MYADML2",  "C2CD4D", "SCARNA7",
                        "EME2", "HIST1H1E", "LGALS9B", "REP15", "RAB19", "BMS1")

#Safe all CRC upregulated genes in a list
crc_upregulated_genes <- c("ABHD8", "ALPL", "BRCC3", "C20orf194", "CANT1", "CASP6",
                           "CCDC106", "CDK8", "CLIP4", "CYB5B", "DFFB", "DFNA5",
                           "DNPEP", "DTX3", "ENGASE", "ETHE1", "EVL", "GNB5", "LBR",
                           "MRAS", "OBSL1", "PCK2", "PPA1", "PREP", "QKI", "RRAGD",
                           "SLC12A2", "SLC26A6", "SMARCA1", "SPG20", "SRPK1", "TCEA2",
                           "TCF7", "UBE2E2", "UCHL3", "ZNF304", "ZNF512B", "ADAT2", "ATP9A",
                           "CALML4", "CNNM4", "DGAT1", "DNAJC15", "FUT4", "GALNT4", "GMDS",
                           "GSTM2", "NDFIP2", "PARP4", "POC1B", "SH3D19", "SNTA1", "TGDS", 
                           "TMEM63A", "WDR91", "KLHDC8B", "MFSD9", "NRARP", "PDZD8", "RNF6",
                           "ROCK2", "STT3B", "ZNF57" )

#Combine all genes in a single list
all_genes <- c(crc_specific_genes, crc_upregulated_genes)

#Select all the cells in the datasets that are epithelial cells
epithelial_seurat <- subset(gse_integrated, subset = Cell_type == "Epithelial cells")

#Make a dotplot of all the genes and select them based on their expression
#If a dot is blue, it has an expression above 0 and if a dot is grey its expression is below 0
#If the dot is big it is expression in more than 25% of the cells and if small less than 25%
DotPlot(epithelial_seurat, features = all_genes, cols=c("light gray", "blue"), col.min = 0, col.max = 0.001, scale.max = 25, scale.min = 24.99)+ 
  coord_flip()+ theme(legend.position = "bottom",legend.box="vertical", axis.title = element_text())+RotatedAxis() + ylab("Cell types") + xlab("Genes")+ scale_size(range = c(1, 3))

#Genes need to show an expression above 0 and be present in 25% of the cells in at least 2/3 tumor epithelial populations to be selected
#Now select the genes that follow this criteria
#te_group contians all the tumor epithelial genes that were selected
te_group <- c("BMS1", "MYADML2", "DUSP27", "TLX1", "LY6G6D", "SLC6A4", "ZNF573",
              "NOXO1", "CEL", "STT3B", "ROCK2", "PDZD8", "TMEM63A", "TGDS", 
              "POC1B", "PARP4", "NDFIP2", "DNAJC15", "DGAT1", "ATP9A", "ZNF512B",
              "UCHL3", "SRPK1", "SLC12A2", "PREP", "PPA1", "DNPEP", "CYB5B",
              "CDK8", "CANT1", "BRCC3")


#================================== Dotplot of tumor epithelial specific genes in all cells ===============================================#

#The next step is to look if any of the genes show expression in any of the stromal cell types
#If they show an expression above 0 in any of the stromal cell types, the genes are not tumor epithelial specific
#if the genes show an expression below 0, they are coloured white and can't be noticed in contrast to the background
DotPlot(gse_integrated, features = te_group, cols=c("white", "blue"), col.min = 0, scale.min = 25) + 
  coord_flip()+ theme(legend.position = "bottom",legend.box="vertical", axis.title = element_blank())+RotatedAxis()+ scale_size(range = c(1, 3))

#There are 13 genes that showed expression in the stromal cell types
leave_out <- c("CDK8", "DNPEP", "ZNF512B", "ATP9A", "DNAJC15", "PARP4", "TGDS", "PDZD8",
               "ROCK2", "STT3B", "ZNF573", "SLC6A4", "BMS1")

#Remove the stromal genes from the te_group
te_group <- te_group[- which(te_group %in% leave_out)]

#Rerun the dotplot with the now fully tumor epithelial specific genes
DotPlot(gse_integrated, features = te_group, cols=c("white", "blue"), col.min = 0, scale.min = 25) + 
  coord_flip()+ theme(legend.position = "bottom",legend.box="vertical", axis.title = element_blank())+RotatedAxis()+ scale_size(range = c(1, 3))
