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
gse132257_seurat <- CreateSeuratObject(counts = gse132257[, epithelial_positions], project = "gse132257", min.cells = 3, min.features = 200)
anno_epi_stro <- anno[epithelial_positions,] #Make an annotation data matrix containing all epithelial and stromal cells
#Add a cancer cell type column in the seurat object (cell disease status + cell type)
gse132257_seurat@meta.data$cancer_cell_type <- paste0(anno_epi_stro$Class,"_",anno_epi_stro$Cell_type)
#Add the cell types to the meta data in the seurat object
gse132257_seurat@meta.data$Cell_type <- anno_epi_stro$Cell_type
#Write the study 
save(gse132257_Seurat,file="gse132257_Seurat.Rda")