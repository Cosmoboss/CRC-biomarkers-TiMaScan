#Import the necessary librbaries
library(mlbench)
library(caret)

#Go to the directory which contains the files
setwd("C:/Users/Gebruiker/Documents/LUMC/secreted_mono_genes")
#Read in all the gene expression data from all cancers
TCGA <- read.table("TCGA_genes_after_filters.xls", sep="\t", header=T)
#Get all of the unique cancers that are listed in the TCGA data
cancers <- unique(TCGA$cancer)


setwd("C:/Users/Gebruiker/Documents/LUMC/TCGA/TCGA_ML")
#Read in all the CRC upregulated genes (with their TCGA name attached)
genes63 <- read.table("63_genes.txt", sep="\t", header=F)

`%notin%` <- Negate(`%in%`) # Create a negate of %in% operator

markers <- c(genes63[c(5, 27, 42, 54),], "cancer")
TCGA_markers <- TCGA[,which(colnames(TCGA) %in% markers)]
#Select all cancers without COADREAD, change them to "Not COADREAD" and
#attach them again to COADREAD
cancer_set <- TCGA_markers[which(TCGA_markers$cancer %notin% c("COADREAD", "COAD", "READ", "ESCA", "STAD")),]
cancer_set$cancer <- "other"
cancer_set <- rbind(cancer_set, TCGA_markers[which(TCGA_markers$cancer %in% c("COADREAD")),])


#Change all "other" to "NO" and all "COADREAD" to "YES"
#I had to do this to be able to run the model.
#It refused to run with the control groups being COADREAD or Not COADREAD
cancer_set[cancer_set == "other"] <- "NO"
cancer_set[cancer_set == "COADREAD"] <- "YES"

cancer_set$cancer <- factor(cancer_set$cancer)

#Make a control training set
control <- trainControl(method = "repeatedcv", number=10, repeats=5, savePredictions = "final", classProbs = T)

#Increase the grid size if it gives an error during RF
#TG = expand.grid(k=1:3,size=seq(5,20,by=5))

#Run RF
model <- train(cancer~., data = cancer_set, method='rf', importance=TRUE, trControl=control)

#Get the variable importance from RF
importance <- varImp(model, scale=FALSE)
plot(importance)

#Make a confusion matrix
Confusion_matrix <- confusionMatrix(model$pred[order(unique(model$pred$rowIndex)),2], cancer_set$cancer)


