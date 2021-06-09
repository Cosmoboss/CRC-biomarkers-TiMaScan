#Load in the necessary libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)

#This code takes the normalized expression values of genes in healthy cell types from the HPA
#The dataset is downloadable at https://www.proteinatlas.org/about/download

#Load in the data from the HPA
HPA_sc_type <- read.table("rna_single_cell_type.tsv", sep="\t", header=T)

#List all the genes that you want to plot in a heatmap
#These are the 18 tumor epithelial specific genes that were selected
te_group <- c("MYADML2", "DUSP27", "TLX1", "LY6G6D", "NOXO1", "CEL", "TMEM63A",
              "POC1B", "NDFIP2", "DGAT1", "UCHL3", "SRPK1", "SLC12A2", "PREP",
              "PPA1", "CYB5B", "CANT1", "BRCC3")

#Select all the rows that include the genes from te_group
HPA_values <- HPA_sc_type[which(HPA_sc_type$Gene.name %in% te_group),]

#Add all unique cell types in a list
Cell_types <- unique(HPA_values$Cell.type)

#Make a matrix using the normalized expression values of all genes and cell types
for(x in 1:length(Cell_types)){
  if(x == 1){
    NX_matrix <- HPA_values[which(HPA_values$Cell.type == Cell_types[x]),]$NX
  }else{
    NX_matrix <- cbind(NX_matrix, HPA_values[which(HPA_values$Cell.type == Cell_types[x]),]$NX)
  } 
}

#Change the column names to the cell types
colnames(NX_matrix) <- Cell_types
#Change the row names to the gene names
row.names(NX_matrix) <- unique(HPA_values$Gene.name)

#The HPA classifies cell types with umbrella term cell types on their website and not in the downloadable data
#By listing all the cell types you can add sub-groups in the heatmap
epithelial_cells <- c("Enterocytes", "Mucus-secreting cells", "Paneth cells",
                      "Ciliated cells", "Club cells", "Exocrine glandular cells",
                      "Basal glandular cells", "Glandular cells", "Basal keratinocytes",
                      "Suprabasal keratinocytes", "Cholangiocytes", "Hepatocytes",
                      "Alveolar cells type 1", "Alveolar cells type 2", "Collecting duct cells",
                      "Distal tubular cells", "Proximal tubular cells", "Ductal cells",
                      "Sertoli cells", "Urothelial cells")

endocrine_cells <- c("Intestinal endocrine cells", "Pancreatic endocrine cells",
                     "Leydig cells")

neuronal_cells <- c("Cone photoreceptor cells", "Rod photoreceptor cells",
                    "Bipolar cells", "Horizontal cells")

glial_cells <- c("Muller glia cells")

germ_cells <- c("Spermatogonia", "Spermatocytes", "Early spermatids", "Late spermatids")

trophoblast_cells <- c("Cytotrophoblasts", "Syncytiotrophoblasts",
                       "Extravillous trophoblasts")

vascular_cells <- c("Endothelial cells")

muscle_cells <- c("Cardiomyocytes", "Smooth muscle cells")

pigment_cells <- c("Melanocytes")

mesenchymal_cells <- c("Fibroblasts", "Ito cells", "Peritubular cells")

undifferentiated_cells <- c("Undifferentiated cells")

blood_and_immune_cells <- c("B-cells", "T-cells", "granulocytes", "Monocytes", "monocytes", "Macrophages",
                            "Hofbauer cells", "Kupffer cells", "Erythroid cells")

#Make a list of grouping
cell_cluster <- colnames(NX_matrix)
cell_cluster[which(cell_cluster %in% epithelial_cells)] <- "epithelial_cells"
cell_cluster[which(cell_cluster %in% endocrine_cells)] <- "endocrine_cells" 
cell_cluster[which(cell_cluster %in% neuronal_cells)] <- "neuronal_cells" 
cell_cluster[which(cell_cluster %in% glial_cells)] <- "glial_cells"
cell_cluster[which(cell_cluster %in% germ_cells)] <- "germ_cells"
cell_cluster[which(cell_cluster %in% trophoblast_cells)] <- "trophoblast_cells"
cell_cluster[which(cell_cluster %in% vascular_cells)] <- "vascular_cells"
cell_cluster[which(cell_cluster %in% muscle_cells)] <- "muscle_cells"
cell_cluster[which(cell_cluster %in% pigment_cells)] <- "pigment_cells"
cell_cluster[which(cell_cluster %in% mesenchymal_cells)] <- "mesenchymal_cells"
cell_cluster[which(cell_cluster %in% undifferentiated_cells)] <- "undifferentiated_cells"
cell_cluster[which(cell_cluster %in% blood_and_immune_cells)] <- "blood_and_immune_cells"

#uncomment the pdf and dev.off() if you want to safe the plot
#pdf("HPA_heatmap_18_genes.pdf", width = 6, height = 8)
Heatmap(scale(t(NX_matrix)),
        split = cell_cluster, #split based on the cell cluster groups
        cluster_row_slices = F,
        row_title_side = "left",
        row_names_side = 'left',
        show_row_dend = TRUE,
        row_names_gp = gpar(fontsize = 7, fontface = 'bold'),
        show_column_dend = TRUE,
        row_title = "Cell types",
        column_title = "Genes",
        column_names_side = "top",
        column_names_rot = 90,
        column_names_gp = gpar(fontsize = 7, fontface = 'bold'),
        name="Gene expression\nper\ncell type\n")
#dev.off()