# Microarray <p1>

In this folder are all the codes that were used to process microarray data.
Do note that the codes do not process all the datasets at once.
However, they are usable for all datasets if you change a few parameters which differ per dataset.

The codes listed are:

* Process_microarray
  * Processes the study GSE62322
* Merge_microarray
  * Merges the datasets from HGU133Plus 2.0 and applies ComBat
* Filtering_Absent_genes_in_microarray
  * Filters out absent/low expressed genes in the microarray platforms by using RNA-Seq genes and Housekeeping genes as reference
