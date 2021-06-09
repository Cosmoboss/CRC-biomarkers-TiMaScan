#import all libraries
library(sva)
library(ggplot2)

#This code uses all studies coming from the microarray platform HGU133Plus 2.0.
#It is applicable for HGU133A as well. Variables however need to be changed however
#Put all the study names in a list. This list is used to go to directories and add study names to samples
studies <- c("gse4107","gse4183", "gse8671", "gse9348", "gse15960", "gse18105", "gse23878", "gse39582", "gse110224")
#Create an empty list with sample amount and study name
gse_meta <- data.frame(samples=character(), gse_study=character())

#Make a variable to count the amount of samples
n_samples = 1

#loop through the list of studies
for (x in 1:length(studies)){
  #Go to the directory with the processed study 
  setwd(paste0("C:/Users/gebruiker/Documents/LUMC/GSE_identifier_script/[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array/", studies[x]))
  #If this is the first instance of the loop, make new variables for the dataframe and pmeta.
  if (x == 1){
    current_df <- read.table(paste0(studies[x], "_topTable.xls"), stringsAsFactors = TRUE)
    current_df = current_df[, 1:(length(current_df)-9)]
    pmeta_all <- read.table(paste0(studies[x], "_pmeta.xls"), header=TRUE)
    merged_df = current_df
  }
  #If this isn't the first instance, combine all dataframes and pmeta's together
  else{
    current_df <- read.table(paste0(studies[x], "_topTable.xls"), stringsAsFactors = TRUE)
    current_df = current_df[, 1:(length(current_df)-9)]
    pmeta_all <- rbind(pmeta_all, read.table(paste0(studies[x], "_pmeta.xls"), header=TRUE))
    merged_df <- cbind(merged_df, current_df)
  }
  #record the total amount of samples and their corrisponding study 
  for (i in 1:length(current_df)){
    meta <- data.frame(samples=n_samples, gse_study=studies[x])
    gse_meta <- rbind(gse_meta, meta)
    n_samples = n_samples+1
  }
  print(studies[x])
}
#Add a columnt to the pmeta with all studies
pmeta_all <- cbind(pmeta_all, study=gse_meta$gse_study)

#Safe the simply merged dataframe
write.table(merged_df, "HGU133Plus2_b4_ComBat_df.xls",sep="\t", quote = F)

#Safe the combined pmeta for the merged dataframe
write.table(pmeta_all, "HGU133Plus2_pmeta.xls",sep="\t", quote = F)
#Do note, currently the dataset is simply merged. ComBat will remove batch effects between studies.
#If you want to see how the samples look simply merged in a PCA plot, skip the combat Steps

#Make a model for combat using the pmeta
model <- model.matrix(~as.factor(ClassMetaData), data=pmeta_all)

#Remove batch effects using ComBat
merged_df <- ComBat(merged_df, batch = pmeta_all$study, mod=model, par.prior=TRUE)

#Safe the merged datasets after batch effects have been removed
write.table(merged_df, "HGU133Plus2_after_Combat_df.xls", sep="\t", quote=F)

#Plotting PCA plot
exprset22=t(merged_df)#invert the table (rows as columns and vice versa)
PC<-prcomp(exprset22) #Run a PCA
PCi<-data.frame(PC$x,Species=pmeta_all$ClassMetaData)#convert the prcomp file to a data frame

#plot the disease types
ggplot(PCi,aes(x=PC1,y=PC2,col=pmeta_all$ClassMetaData))+
  geom_point(size=3)+ #Dot size
  stat_ellipse() + #Confidence interval ellipse at 100%
  stat_ellipse(level=0.80,linetype = 2)+ #Condifence interval ellipse at 80%
  scale_color_manual(values = c("blue", "black"), labels = c("Tumor", "Normal"))+ #Colours and cohorts
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

#plot the studies overlap
ggplot(PCi,aes(x=PC1,y=PC2,col=as.factor(pmeta_all$study)))+
  geom_point(size=3)+
  stat_ellipse() +
  stat_ellipse(level=0.80,linetype = 2)+
  labs(col = "Studies") +
  theme(panel.background = element_blank(),
        axis.line  = element_line(colour = "black"),
        legend.text = element_text(face = "bold"),
        legend.key = element_rect(fill = 'white', size = 0.5, linetype='solid'),
        legend.key.size = unit(1.5, 'lines'),
        #legend.title=element_blank(),
        axis.text.x=element_text(size=rel(1.5), face="bold", colour = "black", angle = 30, hjust=1),
        axis.text.y=element_text(size=rel(1.5), face="bold", colour = "black"),
        axis.title=element_text(size=rel(1.5)))
