# Single-cell <p1>

In this folder are all the codes that were used to process single-cell data and additional files which can be easily loaded into R.
Do note that the codes do not process all the datasets at once.
However, they are usable for all datasets if you change a few parameters which differ per dataset.
The additional files are available on the google drive link below as sometimes R fails to process the files due to the size of the datasets.
  https://drive.google.com/drive/folders/1SKSUQ-e94XoNBqnECdfrOxbvpXfA2WQg?usp=sharing
  Note: The files were to big to upload on GitHub
  
* Make_and_process_seurat_all_cells
  * Reads the count data for the single-cell study and processes it in Seurat for all cell types present
* create_gse132257_epi_stro_seurat
  * Reads in the count data for gse132257 and saves all epithelial and stromal cell types which it processes and writes away in a .Rda file
* Anchor_merge_seurat_objects
  * Merges all seurat objects that contain epithelial and stromal cells together into one big seurat object which can be used for later PCA and dotplots
* Seurat_tumor_epithelial_specific_PCA_&_Dotplots
  * Uses the anchor based merged seurat object to make PCA plots and dotplots which can be used to select tumor epithelial specific genes
* HPA_normalized_expression_heatmap_maker
  * uses the normalized gene expression of 51 healthy cell types in the body and plots the tumor epithelial specific genes against them in a heatmap
