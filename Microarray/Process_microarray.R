#Import all the libraries that will be used
#The cdf and .db libraries depend on the microarray platform 
#use either hgu133a/hgu133plus2 ... hsentrezg ... cdf/.db
#You can install these two packages by downloading the .tar files from brainarray
#and install them via the packages tab in R
library(hgu133ahsentrezgcdf)
library(hgu133ahsentrezg.db)

library(affy)
library(annotate)
library(ggplot2)
library(limma)
library(edgeR)

#Go to the directory where all the files are stored
getwd()
setwd("C:/Users/Gebruiker/Documents/LUMC/GSE_identifier_script/[HG-U133A] Affymetrix Human Genome U133A Array/GSE62322")


#Save all .CEL files from the current directory in a variable called gse_files
gse_files <- list.celfiles()
#Make file names using the GSE identifier in variables for later use
identifier = "gse62322"
identifier_raw = paste0(identifier, "-raw")
identifier_voom = paste0(identifier, "-voom")
identifier_processed_xls = paste0(identifier, "_processed.xls")

#Make an affybatch using the .CEL files and corrisponding Content Definition File (cdf)
#The cdfname can either be HGU133Plus2_Hs_ENTREZG or HGU133A_Hs_ENTREZG
gse_aBatch <- ReadAffy(filenames = gse_files, cdfname = "HGU133A_Hs_ENTREZG")
gse.eset.mas=mas5(gse_aBatch) #Use MAS5 to correct the probes
gse.eset.mas=exprs(gse.eset.mas) #Make an expression set.
gse_safe <- gse.eset.mas #Safe the expression set in an additional variable for later use

gse.eset.mas <- gse_safe[,- c(20,31)] #Remove potential outliers based on sample position
write.table(gse.eset.mas, identifier_raw, sep="\t") #Save the current expression matrix in a file

#Make a pmeta data matrix corrisponding to the sample files
#Do note that this code for the pmeta relies on the tumor and normal samples being
#clustered together in the directory based on tissue and not randomized.
#Two of the same loops to make the C/N list are commented. 
#The difference is that the "N" and "C" are switched.
#You can comment/uncomment the code depending on which position the samples are listed in the files

c <- 21 #amount of cancer tissue samples
n <- 17 #amount of normal tissue samples

samples_list <- c() #Make a list with sample tissues and numbers

#If your files are formatted with starting with normal and afterwards cancer use the two for loops below

#for (x in 1:n){samples_list[x] <- paste0("N",x)}
#for (x in 1:c){samples_list[length(samples_list)+1] <- paste0("C",x)}

#If your files are formatted with starting with cancer and afterwards normal use the two for loops below

for (x in 1:c){samples_list[x] <- paste0("C",x)}
for (x in 1:n){samples_list[length(samples_list)+1] <- paste0("N",x)}

ClassMetaData_list <- c() #Make a list contained either normal/cancer for each sample

#If your files are formatted with starting with normal and afterwards cancer use the two for loops below

#for (x in 1:n){ClassMetaData_list[x] <- "normal"}
#for (x in 1:c){ClassMetaData_list[length(ClassMetaData_list)+1] <- "cancer"}

#If your files are formatted with starting with cancer and afterwards normal use the two for loops below

for (x in 1:c){ClassMetaData_list[x] <- "cancer"}
for (x in 1:n){ClassMetaData_list[length(ClassMetaData_list)+1] <- "normal"}

#Construct a metadata file with sample name, sample tissue and its corrisponding data file
pmeta_safe <- data.frame(samples=samples_list, ClassMetaData=ClassMetaData_list, GSM_Files=gse_files)
#Save it once for later potential comparisons for before/after outlier removal.
write.table(pmeta_safe, file=paste0(identifier,"_pmeta_before_outliers.xls"), sep = "\t", quote=FALSE, row.names = FALSE)

pmeta <- pmeta_safe[- c(20,31),] #Remove any potential outliers if needed
#or
pmeta <- pmeta_safe #If there are no outliers (yet) , use this

#Safe the pmeta file
write.table(pmeta, file=paste0(identifier,"_pmeta.xls"), sep = "\t", quote=FALSE, row.names = FALSE)




#Read in the expression set
gse_raw <- read.table(identifier_raw, sep="\t", header=TRUE, stringsAsFactors = FALSE )
mod<-model.matrix(~0+as.factor(pmeta$ClassMetaData)) #Make a model 
#log2 transform
eset=log2(gse_raw)
# Normalize the raw matrix using vooma {limma}
eset.mas.voom=vooma(eset,design = mod, plot = TRUE)

#Safe the log2 and normalized expression values in a table
write.table(eset.mas.voom$E, identifier_voom, sep="\t")
#Safe the same values in a variable
voom = eset.mas.voom$E

#Plotting PCA plot
#If there are outliers present in the data you can manually pick them here and remove them
#If an outlier is removed, start from code line 33 until here again while removing the outlier(s) 
exprset22=t(voom) #invert the table (rows as columns and vice versa)
PC<-prcomp(exprset22) #Run a PCA plot
PCi<-data.frame(PC$x,Species=pmeta$ClassMetaData) #convert the prcomp file to a data frame
#pdf("GSE62322-pca.pdf", width=5, height=5) #if you want to save a pdf, run this and the ggplot code below and finish with dev.off()
ggplot(PCi,aes(x=PC1,y=PC2,col=as.factor(pmeta$ClassMetaData)))+ #Load in the data
  geom_point(size=3)+ #Dot size
  geom_text(size=5, aes(label=pmeta$samples),hjust = 0, nudge_x = 0.1) + #Text size
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
#dev.off()

#Adding gene names to the table using the corrisponding database
Gene <- getSYMBOL(rownames(voom), "hgu133ahsentrezg.db")
results <- cbind(voom, Gene) #add a gene column to the result

#Safe the processed microarray file with corrisponding gene names.
write.table(results, file = identifier_processed_xls, sep = "\t", quote=FALSE)
