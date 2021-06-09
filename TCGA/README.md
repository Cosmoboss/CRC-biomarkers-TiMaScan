# TCGA <p1>

In this folder are all the codes that were used to process, filter and merge the raw cancer data from TCGA downloaded on GDAC Broad Firehose.
In addition there are codes for creating a correlation plot, PCA plot and run a random forest machine learning model using the data.

The codes listed are:

* TCGA_read_filter_and_normalize
  * Processes and filters the raw TCGA data
* TCGA_merge_&_pmeta_&_filter_secr_macr
  * merges all the datasets, constructs a huge pmeta and filters the genes for secreted proteins and macrophage/monocyte specific genes
* TCGA_topTable_maker_&_upregulated_genes
  * Does differential expression between all other cancer types and makes an presence/absence matrix to look for genes that are upregulated in CRC 
* TCGA_selected_markers_PCA
  * Makes a PCA plot containing the expression values from CRC upregulated genes in all cancer types to see how they differentiate CRC from other cancer types
* Markers_Random_Forest
  * Runs a Random Forest model to see if the final set of markers differentiate CRC samples from all other cancer types
* Genes_in_cancers_means_correlation_plot
  * Makes a correlation plot of average gene expression per cancer type to see how other cancers can show similar expression to CRC
