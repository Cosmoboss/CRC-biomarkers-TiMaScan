#Load in the libraries
library(Seurat)

#Load in the the count data from the single-cell study
gse132257 <- read.table("GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt", sep="\t", header = TRUE)
row.names(gse132257) <- gse132257$Index #Make the first column (gene names) the row names
gse132257 <- gse132257[,- 1] #Remove the first column
#Load in the annotation file from the single-cell study
anno <- read.table("GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt", sep="\t", header = TRUE)
#Get all positions from the annotation table that are epithelial or stromal cells
epithelial_positions <- as.numeric(row.names(anno[which(anno$Cell_type %in% c("Epithelial cells", "Stromal cells")),]))

#The positions for the samples in the annotation file are in the same order as the column names for the count data.
#Make a seurat object with the epithelial and stromal cells from the count data
gse132257_seurat <- CreateSeuratObject(counts = gse132257, project = "gse132257", min.cells = 3, min.features = 200)
#Add a cancer cell type column in the seurat object (cell disease status + cell type)
gse132257_seurat@meta.data$cancer_cell_type <- paste0(anno$Class,"_",anno$Cell_type)
#Add the cell types to the meta data in the seurat object
gse132257_seurat@meta.data$Cell_type <- anno_epi_stro$Cell_type

#Calculate the percentage of mitochondiral genes in the datasets
gse132257_seurat[["percent.mt"]] <- PercentageFeatureSet(gse132257_seurat, pattern = "^MT-")
#a high percentage of mitochondrial genes suggests the cell is of bad quality.
#Remove any gene with more than 5% of mitochondiral genes
gse132257_seurat <- subset(gse132257_seurat, subset = percent.mt < 5)
#Normalize the count table in seurat
gse132257_seurat <- NormalizeData(gse132257_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#Find the top 2000 variable features (genes) in the seurat object
gse132257_seurat <- FindVariableFeatures(gse132257_seurat, selection.method = "vst", nfeatures = 2000)
#Scale the data
gse132257_seurat <- ScaleData(gse132257_seurat, verbose = FALSE)
#Plot a PCA
gse132257_seurat <- RunPCA(gse132257_seurat, npcs = 30, verbose = FALSE)
#Make a UMAP
gse132257_seurat <- RunUMAP(gse132257_seurat, reduction = "pca", dims = 1:30)
#Plot the UMAP
DimPlot(gse132257_seurat, reduction = "umap", group.by = "Cell_type")
