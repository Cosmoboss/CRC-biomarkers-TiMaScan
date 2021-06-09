#Load the libraries that will be used.
library(matrixStats)

#Read in the merged microarray datasets
hgu133plus2 <- read.table("HGU133Plus2_after_ComBat_df.xls", sep="\t", header=T)
hgu133a <- read.table("HGU133A_after_ComBat_df.xls", sep="\t", header=T)

#Read in the RNA-Seq studies
gse74369 <- read.table("gse74369_processed.xls", header=T)
gse146009 <- read.table("gse146009_processed.xls", sep = "\t", header=T)

#Remove any gene that is an NA
gse74369 <- gse74369[!is.na(gse74369$external_gene_name),]
gse146009 <- gse146009[!is.na(gse146009$external_gene_name),]

#Identify all unique genes in RNA-Seq
rna_seq_genes <- unique(c(gse74369$external_gene_name, gse146009$external_gene_name))

`%notin%` <- Negate(`%in%`) # Create a negate of %in% operator

hgu133plus2_genes <- unique(hgu133plus2$Gene) #All unique genes in HGU133Plus2
hgu133a_genes <- unique(hgu133a$Gene) #All unique genes in HGU133A

#Find all HGU133's genes not present in RNA-Seq genes
hgu133a_specific <- hgu133a_genes[hgu133a_genes %notin% rna_seq_genes]
hgu133plus2_specific <- hgu133plus2_genes[hgu133plus2_genes %notin% rna_seq_genes]

#Overlap both HGU133 to find microarray specific genes
microarray_specific <- hgu133a_specific[hgu133a_specific %in% hgu133plus2_specific]
#find all genes specific in HGU133plus2 only
hgu133plus2_specific <- hgu133plus2_specific[hgu133plus2_specific %notin% microarray_specific]

#List of the Housekeeping genes
HKG_genes <- c("GPI","PSMB2","PSMB4","SNRPD3","RAB7A","REEP5","VCP","CHMP2A","C1orf43","VPS29","EMC7","CSTB","GUSB","HADHA","HPRT1",
               "SGSH","PGK1","TCOF1","CYB5R3","GM2A","SOD1","GNAS","COMT","GUK1","IMPDH2","P4HB","POLR2A","RPL5","RPL8","RPL11","RPL19",
               "RPL17","RPL27","RPL32","RPL34","RPL36AL","RPS5","RPS9","RPS11","RPS13","RPS15","RPS16","RPS24","RPS25","ADAR","ADD1",
               "AP1B1","AES","ANXA6","ATP6AP1","BTF3","ENTPD6","GAPDH","CHD4","TBCB","CSNK2B","CTBP1","DAD1","DAXX","DDT","EIF4G2",
               "ENO1","FBL","EXTL3","XRCC6","BLOC1S1","GDI1","GDI2","PRMT1","HSBP1","APLP2","ARAF","ARF1","ARF4","ARF5","RHOA","ATF4",
               "ATP5D","ATP5G3","ATP6V0C","ATP6V1E1","ATP5O","BSG","CANX","CAPNS1","SEPT7","CENPB","CLTA","CLTB","COX4I1","COX5B",
               "COX6B1","COX7A2","COX7C","CTNNB1","CTSD","CYC1","E2F4","EEF2","FAU","FTH1","GAPDH","GOT2","GPX4","HMGB1","HNRNPD",
               "HNRNPK","ID3","JAK1","M6PR","MAP4","MAZ","MIF","MAP3K11","MPG","MTX1","NDUFA2","NDUFC1","NME2","ODC1","PRDX1","PFDN1",
               "PFDN5","SLC25A3","PHF1","PIM1","PPP1R10","PRKAG1","PRKCSH","PSMA7","PSMB1","PSMB7","PSMD2","PSMD3","PSMD8","PSMD11",
               "PSME2","PTBP1","RING1","RNH1","RPA2","RPL15","RPN1","RPS27A","SAFB","SRSF2","SGTA","SNRNP70","SNRPB","SNRPG","SRM",
               "SRP14","SSR2","TAPBP","TMBIM6","TSTA3","TTC1","TUFM","TXN","UBA1","UBE2D2","UBE2I","UQCRC1","YWHAB","YWHAZ","ZNF91",
               "HIST1H2BC","SLC25A11","STK24","YARS","EIF3D","EIF3F","EIF3G","EIF3I","BECN1","B4GALT3","SNX3","GPAA1","ADAM15","BANF1",
               "ARHGEF7","MCM3AP","BUD31","CPNE1","RPS6KB2","UBE2M","RPL14","ATP5A1","ATP6V0B","AP2M1","AP2S1","COX8A","NDUFB7","RAB1A",
               "SDHA","SREBF1","ATP6V1F","ARHGAP1","ARHGDIA","PTTG1IP","CD81","SEPT2","ENSA","ERH","HDGF","HNRNPAB","ILF2","ILK",
               "NDUFA1","NDUFS5","SNRPA","SNRPD2","DDX39B","PABPN1","C21orf33","MTA1","HGS","COX7A2L","MAPKAPK2","VAMP3","ATP6V1G1",
               "ATP5J2","SPAG7","H2AFY","C14orf2","CLOCK","ACTN4","CAPZB","FUS","NDUFA7","PFN1","RBM8A","WDR1","ABL1","ATP5G1",
               "BMI1","DDOST","GNB2","HINT1","HSPA5","JUND","NCL","MLF2","CFL1","HNRNPH1","KARS","LAMP1","LDHA","RPS14","SCAMP3",
               "ARPC4","ARPC3","TSFM","ARPC2","BCAP31","TRIM28","EIF1","SRRM1","SAP18","ACAT2","SLC25A1","VPS72","CCT3","UQCRFS1",
               "UQCRH","RPL10","AKR1A1","TUBA1B","HAX1","NEDD8","PIN1","DPF2","VARS","RAN","C1D","ZNHIT1","TIMM44","TADA3","ATP5H",
               "NXF1","CREB3","SYNCRIP","HYOU1","ANP32B","AGPAT1","RABAC1","CCT7","DRAP1","PRPF8","HEXIM1","TRIM27","SARS","RRAGA",
               "API5","HSPA8","TUBGCP2","JTB","NUDT3","RNPS1","TALDO1","ZFPL1","AFG3L2","KDELR1","SEC61B","TMED2","YWHAQ","UQCR11",
               "COPS6","CALM1","IDH3B","RAC1","SUMO3","RTN4","KAT7","ATP5I","NDUFV1","RPL10A","TCEB2","RPL35","ATXN2L","LYPLA2",
               "PARK7","COPE","GABARAP","GABARAPL2","ABL1","HSP90AB1","CASC3","NONO","CD3EAP","DNPEP","ARL2BP","AHSA1","CIZ1","AATF",
               "FBXO7","PICK1","H2AFV","RPL13A","PDCD6","EIF3K","PRRC2B","PPP2R1A","CNPY2","PUF60","SEC61G","SND1","UQCRQ","ZNF592",
               "MLEC","PTDSS1","IST1","EFCAB14","MFN2","PDAP1","LMTK2","TCF25","XPO7","ESYT1","CTDNEP1","BRMS1","NELFB","RAP1B","ANAPC5",
               "TRAP1","INPP5K","PTOV1","TMED9","OTUB1","BTBD2","COMMD4","UBB","TERF2IP","TOMM7","GSK3A","SAR1A","STARD7","SDR39U1",
               "UBC","NDUFV2","MRPS12","POLR2L","MRPL23","PPP1R11","MCL1","POLR2F","RELA","TUT1","CDK11A","KXD1","TMEM109","C9orf16",
               "MAP2K2","MRPL9","TMEM147","MYL12B","ZNF384","TEX261")

#put in the microarray platform that you want to filter
a1 <- hgu133plus2

#Remove any genes that are labelled as NA
a1 <- a1[!is.na(a1$Gene),]

#Calculate the average expression of the genes in all samples
a1$AveExpr <- rowMeans(a1[,1:(length(a1)-1)])

#Add this if you're curious about the hgu133plus2 specific genes
sel_3381 <- a1[a1$Gene %in% hgu133plus2_specific,]

#Safe the expression levels of the microarray specific genes
sel_1504 <- a1[a1$Gene %in% microarray_specific,]

columns_sel <- c("Gene","AveExpr") #Save these columns for later

#get the selected columns from hgu133plus2 specific genes if you're processing it
sel_3381_sel <- sel_3381[,columns_sel]

#Safe the gene names and average expression of the microarray specific genes
sel_1504_sel <- sel_1504[,columns_sel]

#Get all HKG genes from all the genes
sel_hkg <- a1 [a1$Gene %in% HKG_genes,]

#Safe the gene names and average expression from the HKG
sel_hkg_sel <- sel_hkg[,columns_sel]

#Plot the average expression of all genes
hist(a1$AveExpr, col="skyblue", ylim = c(0, 2500),border=F) #All genes in microarray
hist(sel_3381$AveExpr, col="gold",add=TRUE,border=F) #hgu133plus2 specific genes
hist(sel_1504$AveExpr, col="red",add=TRUE,border=F) #microarray specific genes
hist(sel_hkg$AveExpr, col="black",add=TRUE,border=F) #Housekeeping genes

#To calculate Standard Deviation (SD)
aa1 <- transform(as.matrix(a1[, 1:(length(a1)-2)]), SD=apply(as.matrix(a1[, 1:(length(a1)-2)]),1, sd, na.rm = TRUE))
A1_Sel <- a1[,columns_sel]
A1_Sel$SD <- aa1$SD

#Add the SD for each GENE type
sel_3381_sel$SD <- A1_Sel[which(row.names(sel_3381_sel) %in% row.names(A1_Sel)),]$SD
sel_1504_sel$SD <- A1_Sel[which(row.names(sel_1504_sel) %in% row.names(A1_Sel)),]$SD
sel_hkg_sel$SD <- A1_Sel[which(row.names(sel_hkg_sel) %in% row.names(A1_Sel)),]$SD

#Calculate the average SD and average expression to identify the cut-off
AveExpr <- mean(apply(sel_1504[,1:(length(sel_1504)-2)], 1, mean))
AveSD <- mean(apply(sel_1504[,1:(length(sel_1504)-2)], 1, sd))

#Test the filtering criteria of AveExpr + SD times 1,2, or 3
#SD less than 1
FilteredGenes_V1 <- A1_Sel[which(A1_Sel$AveExpr < AveExpr + (AveSD * 1)),]
SelectedGenes_V1 <- A1_Sel[- which(A1_Sel$AveExpr < AveExpr + (AveSD * 1)),]
filtered_1504_V1 <- sel_1504_sel[which(row.names(sel_1504_sel) %in% row.names(FilteredGenes_V1)),]
filtered_HKG_V1 <- sel_hkg_sel[which(row.names(sel_hkg_sel) %in% row.names(FilteredGenes_V1)),]
print(c(dim(filtered_1504_V1), dim(filtered_HKG_V1), dim(SelectedGenes_V1), (AveExpr+AveSD*1)))

#SD less than 2
FilteredGenes_V2 <- A1_Sel[which(A1_Sel$AveExpr < AveExpr + (AveSD * 2)),]
SelectedGenes_V2 <- A1_Sel[- which(A1_Sel$AveExpr < AveExpr + (AveSD * 2)),]
filtered_1504_V2 <- sel_1504_sel[which(row.names(sel_1504_sel) %in% row.names(FilteredGenes_V2)),]
filtered_HKG_V2 <- sel_hkg_sel[which(row.names(sel_hkg_sel) %in% row.names(FilteredGenes_V2)),]
print(c(dim(filtered_1504_V2), dim(filtered_HKG_V2), dim(SelectedGenes_V2), (AveExpr+AveSD*2)))

#SD less than 3
FilteredGenes_V3 <- A1_Sel[which(A1_Sel$AveExpr < AveExpr + (AveSD * 3)),]
SelectedGenes_V3 <- A1_Sel[- which(A1_Sel$AveExpr < AveExpr + (AveSD * 3)),]
filtered_1504_V3 <- sel_1504_sel[which(row.names(sel_1504_sel) %in% row.names(FilteredGenes_V3)),]
filtered_HKG_V3 <- sel_hkg_sel[which(row.names(sel_hkg_sel) %in% row.names(FilteredGenes_V3)),]
print(c(dim(filtered_1504_V3), dim(filtered_HKG_V3), dim(SelectedGenes_V3), (AveExpr+AveSD*3)))


#Safe the criteria where the least amount of HKG were filtered whilst filtering
#as much absent/low expressed genes (SD times 1 in my case)
safe_filtered <- a1[which(row.names(a1) %in% row.names(FilteredGenes_V1)),]
write.table(safe_filtered, "HGU133Plus2_filtered.xls", sep="\t", quote=F)
safe_selected <- a1[which(row.names(a1) %in% row.names(SelectedGenes_V1)),]
write.table(safe_selected, "HGU133Plus2_after_lowExpression_filter.xls", sep="\t", quote=F)
