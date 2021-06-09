#Import all the packages that will be used
library(edgeR)
library(sva)
library(ggplot2)

#Go to the directory with the processed pmeta for each study
setwd("C:/Users/Gebruiker/Documents/LUMC/GSE_identifier_script/RNA-Seq/gse74369")
p1 <- read.table("gse74369_pmeta_after_outliers.xls", sep="\t", header = TRUE)

setwd("C:/Users/Gebruiker/Documents/LUMC/GSE_identifier_script/RNA-Seq/gse146009")
p2 <- read.table("gse146009_pmeta_after_outliers.xls", sep="\t", header = TRUE)

#Add another column corrisponding to the study and merge the pmetas
p1$Study <- "gse74369"
p2$Study <- "gse146009"
pmeta_all <- rbind(p1,p2)

#Change the rownames to a logical order
rownames(pmeta_all) <- 1:nrow(pmeta_all)

#Go to a directory with all the count files for all studies
setwd("C:/Users/Gebruiker/Documents/LUMC/GSE_identifier_script/RNA-Seq/RNA_seq_countfiles")
#Read in all the counts by using the sample names from the pmeta
counts1 = readDGE(pmeta_all$countf)$counts
noint = rownames(counts1) %in%
  c("__no_feature","__ambiguous","__too_low_aQual",
    "__not_aligned","__alignment_not_unique") #Safe these rownames to remove them later
cpms = cpm(counts1) #Calculate the CPMS
keep = rowSums(cpms>1)>=5 & !noint #If a sample has a CPM above 1 in 5 or more samples it isn't absent/low expressed
counts = counts1[keep,] #Keep all the counts that were present/high abundant

#Run Combat-Seq
adjusted <- ComBat_seq(counts, batch=as.factor(pmeta_all$Study), group=pmeta_all$ClassMetaData)

#Run PCA
exprset22=t(adjusted)
PC<-prcomp(exprset22)
PCi<-data.frame(PC$x,Species=pmeta_all$ClassMetaData)

#PCA for the disease types
ggplot(PCi,aes(x=PC1,y=PC2,col=pmeta_all$ClassMetaData))+
  geom_point(size=3)+ #Dot size
  #geom_text(size=5, aes(label=pmeta$samples),hjust = 0, nudge_x = 0.1) + #Text size
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

#PCA for the studies
ggplot(PCi,aes(x=PC1,y=PC2,col=as.factor(pmeta_all$Study)))+
  geom_point(size=3)+ #Dot size
#  geom_text(size=5, aes(label=pmeta$samples),hjust = 0, nudge_x = 0.1) + #Text size
  stat_ellipse() + #Confidence interval ellipse at 100%
  stat_ellipse(level=0.80,linetype = 2)+ #Condifence interval ellipse at 80%
  scale_color_manual(values = c("gold", "forest green"), labels = c("GSE146009", "GSE74369"))+ #Colours and cohorts
  labs(col = "Study") + #Legend title
  theme(panel.background = element_blank(), #Extra bolding and positioning of x/y-axis
        axis.line  = element_line(colour = "black"),
        legend.text = element_text(face = "bold"),
        legend.key = element_rect(fill = 'white', size = 0.5, linetype='solid'),
        legend.key.size = unit(1.5, 'lines'),
        #legend.title=element_blank(),
        axis.text.x=element_text(size=rel(1.5), face="bold", colour = "black", angle = 30, hjust=1),
        axis.text.y=element_text(size=rel(1.5), face="bold", colour = "black"),
        axis.title=element_text(size=rel(1.5)))
