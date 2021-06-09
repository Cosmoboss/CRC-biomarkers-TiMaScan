#Import the necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

`%notin%` <- Negate(`%in%`) # Create a negate of %in% operator

#Load the seurat objects that include only the epithelial and stromal cells
load("~/GSE132257/gse132257_seurat.Rda")
load("~/GSE132465/gse132465_seurat.Rda")
load("~/GSE144735/gse144735_seurat.Rda")

#We want cell types that are common in all three single-cell studies
#Select only the Normal and Tumor cell types from gse144735 and exclude the border cell types
gse144735_seurat_sub <- subset(x = gse144735_seurat, subset = cancer_cell_type %notin% c("Border_Epithelial cells", "Border_Stromal cells"))

#Calculate the amount of mitochondrial genes there are present in the datasets
gse132257_seurat[["percent.mt"]] <- PercentageFeatureSet(gse132257_seurat, pattern = "^MT-")
gse132465_seurat[["percent.mt"]] <- PercentageFeatureSet(gse132465_seurat, pattern = "^MT-")
gse144735_seurat_sub[["percent.mt"]] <- PercentageFeatureSet(gse144735_seurat_sub, pattern = "^MT-")

#A high percentage of mitochondrial genes suggests the cell is dead or of bad quality
#Therefore, remove any cell with more than 5% mitochondrial genes
gse132257_seurat <- subset(gse132257_seurat, subset = percent.mt <5)
gse132465_seurat <- subset(gse132465_seurat, subset = percent.mt < 5)
gse144735_seurat_sub <- subset(gse144735_seurat_sub, subset = percent.mt < 5)

#Make a list with all seurat objects and add their corrisponding study names to the list
gse_list <- list(gse132257_seurat, gse132465_seurat, gse144735_seurat_sub)
names(gse_list) <- c("gse132257", "gse132465", "gse144735")

#Loop through list normalize to normalize the seurat object and find the top 2k varaible features for each seurat object
for (i in 1:length(gse_list)) {
  gse_list[[i]] <- NormalizeData(gse_list[[i]], verbose = FALSE)
  gse_list[[i]] <- FindVariableFeatures(gse_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

#Find a set of anchors between the list of Seurat objects. 
#These anchors can later be used to integrate the objects using the IntegrateData function.
gse_anchors <- FindIntegrationAnchors(object.list = gse_list, dims = 1:30)

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

#Perform a dataset integration using the pre-computed anchors in gse_anchors 
#Integrate all anchors and also all possible CRC genes that weren't included already
gse_integrated <- IntegrateData(anchorset = gse_anchors, dims = 1:30, features.to.integrate = unique(c(gse_anchors@anchor.features, all_genes)))

#Change the default assay to the integrated counts instead of the raw 
DefaultAssay(gse_integrated) <- "integrated"


#Scale and center features from the datasets
gse_integrated <- ScaleData(gse_integrated, verbose = FALSE)
#Run a PCA dimensionality reduction on the seurat object
#This is needed to run a UMAP
gse_integrated <- RunPCA(gse_integrated, npcs = 30, verbose = FALSE)
#Run a UMAP dimensional reduction technique on the seurat object
gse_integrated <- RunUMAP(gse_integrated, reduction = "pca", dims = 1:30)
#Plot all cells from the seurat object in an UMAP whilst grouped by cancer cell type
DimPlot(gse_integrated, reduction = "umap", group.by = "cancer_cell_type")


#Now we add a short code to the merged dataset for each their corrisponding study
#Put all unique orig.idents in a list
orig.idents <- unique(gse_integrated@meta.data$orig.ident)
FRZ <- orig.idents[1:10] #GSE132257
SMC <- orig.idents[11:43] #GSE132465
KUL <- orig.idents[44:55] #GSE144735

#change the orig.idents with the 3-letter tags representing  each study 
gse_integrated@meta.data[which(gse_integrated@meta.data$orig.ident %in% FRZ),1] <- "FRZ"
gse_integrated@meta.data[which(gse_integrated@meta.data$orig.ident %in% SMC),1] <- "SMC"
gse_integrated@meta.data[which(gse_integrated@meta.data$orig.ident %in% KUL),1] <- "KUL"

#Make a new column in the data representing studies with cell types
gse_integrated[["study_cell_type"]] <- paste(gse_integrated@meta.data[,1], gse_integrated@meta.data[,4])

#when you run Idents(gse_integrated) all the levels are still corrisponding to the orig.idents even though they have been removed
#This will change the levels to the study_cell_type row
Idents(gse_integrated) <- factor(gse_integrated$study_cell_type)

#Now that the study has been processed, write it away in a seperate file for later use
save(gse_integrated, file="gse_integrated.Rda")
