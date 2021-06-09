library(ggplot2)
library(limma)


#Read in all the gene expression data from all cancers
TCGA <- read.table("TCGA_genes_after_filters.xls", sep="\t", header=T)
#Get all of the unique cancers that are listed in the TCGA data
cancers <- unique(TCGA$cancer)

#Read in all the CRC upregulated genes (with their TCGA name attached)
genes_TCGA <- read.table("TCGA_63_genes.txt", sep="\t", header=F)

`%notin%` <- Negate(`%in%`) # Create a negate of %in% operator

#Make a cancer set wherein CRC, ESCA and STAD are listed and all other cancers as "other"
cancer_set <- TCGA[which(TCGA$cancer %notin% c("COADREAD", "COAD", "READ", "ESCA", "STAD")),]
cancer_set$cancer <- "other"
cancer_set <- rbind(cancer_set, TCGA[which(TCGA$cancer %in% c("COADREAD", "ESCA", "STAD")),])


#The PCA plot showing how ESCA/STAD are similar to CRC
exprset22=cancer_set[,which(colnames(cancer_set) %in% genes_TCGA$V1)]
PC<-prcomp(exprset22)
PCi<-data.frame(PC$x,Species=factor(cancer_set$cancer))

#pdf("test_PCA.pdf", width=7, height=5) #if you want to save the PCA plot
ggplot(PCi,aes(x=PC1,y=PC2,col=factor(cancer_set$cancer)))+
  geom_point(size=1, alpha=2)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue", "black", "grey", "red"))+
  labs(col = "Cohort") +
  stat_ellipse() +
  stat_ellipse(level=0.80,linetype = 2)+
  theme(panel.background = element_blank(),
        axis.line  = element_line(colour = "black"),
        legend.text = element_text(face = "bold"),
        legend.key = element_rect(fill = 'white', size = 0.5, linetype='solid'),
        legend.key.size = unit(1.5, 'lines'),
        #legend.title=element_text("TCGA cancers comparison"),
        axis.text.x=element_text(size=rel(1.5), face="bold", colour = "black", angle = 30, hjust=1),
        axis.text.y=element_text(size=rel(1.5), face="bold", colour = "black"),
        axis.title=element_text(size=rel(1.5)))
#dev.off()

#Make a cancer set with CRC and all other cancers excluding ESCA/STAD
cancer_set <- TCGA[which(TCGA$cancer %notin% c("COADREAD", "COAD", "READ", "ESCA", "STAD")),]
cancer_set$cancer <- "other"
cancer_set <- rbind(cancer_set, TCGA[which(TCGA$cancer %in% c("COADREAD")),])


#The PCA plot with the selected markers and how they differently classify CRC as compared to other cancer
exprset22=cancer_set[,which(colnames(cancer_set) %in% c("CANT1.124583", "SLC12A2.6558", "DGAT1.8694", "DGAT1.8694"))]
PC<-prcomp(exprset22)
PCi<-data.frame(PC$x,Species=factor(cancer_set$cancer))

ggplot(PCi,aes(x=PC1,y=PC2,col=factor(cancer_set$cancer)))+
  geom_point(size=1, alpha=2)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue",  "grey"))+
  labs(col = "Cohort") +
  stat_ellipse() +
  stat_ellipse(level=0.80,linetype = 2)+
  theme(panel.background = element_blank(),
        axis.line  = element_line(colour = "black"),
        legend.text = element_text(face = "bold"),
        legend.key = element_rect(fill = 'white', size = 0.5, linetype='solid'),
        legend.key.size = unit(1.5, 'lines'),
        #legend.title=element_text("TCGA cancers comparison"),
        axis.text.x=element_text(size=rel(1.5), face="bold", colour = "black", angle = 30, hjust=1),
        axis.text.y=element_text(size=rel(1.5), face="bold", colour = "black"),
        axis.title=element_text(size=rel(1.5)))
