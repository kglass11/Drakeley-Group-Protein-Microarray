#Sanger Data Analysis Script
#Katie Glass
#updated: 9/11/18

#####################################
############DATA ANALYSIS############
#####################################

rm(list=ls())

#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens/Sanger NM V2"
#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens"
#
setwd("/Users/Katie/Desktop/R files from work/100817 Sanger/Sanger NM V2")
getwd()

require("gtools")

library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(reshape2)
library(corrplot)

load(file="Sanger.2.Update.RData")
load(file = "sangerNMcutoffsfinal.RData")

###### Isolate data for non-malarial antigens and test samples only (no controls) or test samples and controls
    ###negative values are included here
    ###this first section is repeated from FMM script (Sanger Analysis Non-Malarial Antigens)

#Remove control protein targets
#Remove samples that should be excluded and isolate test samples - leaves 1325 samples
#This means that the two neg pools and 6 malaria positive controls are not affecting the FMMs
subdata <- norm3.matrix[-rmsamp_all,!(colnames(norm3.matrix) %in% samples_exclude)]
testdata <- subdata[,(colnames(subdata) %in% samples_test)]

#Replace current target names with original target names now that control targets are removed
testdata.df <- merge(testdata, annotation_targets.df, by ="row.names", sort = FALSE)
testdata.df <- tibble::column_to_rownames(testdata.df, var="Row.names")
row.names(testdata.df) <- testdata.df$Name
Ftestdata <- testdata.df[,1:ncol(testdata)]

#do the same thing for the matrix that has the controls - need this later
#Replace current target names with original target names now that control targets are removed
subdata.df <- merge(subdata, annotation_targets.df, by ="row.names", sort = FALSE)
subdata.df <- tibble::column_to_rownames(subdata.df, var="Row.names")
row.names(subdata.df) <- subdata.df$Name
Fsubdata <- subdata.df[,1:ncol(subdata)]

#Merge with target metadata to filter based on expression tag etc.
target.df <- merge(target_meta.df, Ftestdata, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
targetsub.df <- merge(target_meta.df, Fsubdata, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

#isolate data for non-malarial antigens
NMdatameta <- filter(target.df, Category == "non_malarial")
NMdata <- NMdatameta[,(ncol(target_meta.df)+1):ncol(NMdatameta)]
rownames(NMdata) <- NMdatameta$Name

min(NMdata) #-7.300968
max(NMdata) #8.524411

#repeat for subdata (test and control samples)
NMallsampmeta <- filter(targetsub.df, Category == "non_malarial")
NMallsampdata <- NMallsampmeta[,(ncol(target_meta.df)+1):ncol(NMallsampmeta)]
rownames(NMallsampdata) <- NMallsampmeta$Name

#Transpose data so antigens are columns - it's now a matrix
tNMdata <- t(NMdata)

tNMallsampdata <- t(NMallsampdata)

##### Remove antigens which we are not analyzing from cutoff matrix and data matrix
rmant <- c("CT706","Pertussis JNIH-5 [100] *", "CT110", "Pertussis JNIH-3 [1] *", "Pertussis JNIH-3 [0.1] *")

finalcut <- as.matrix(cutoffsavedfinal[!rownames(cutoffsavedfinal) %in% rmant,])

itta <- as.data.frame(tNMdata[,!colnames(tNMdata) %in% rmant])

############### Correlation between antigens #################

### do correlation plots with ALL data, including negative values. 

#run correlation significance test with pearson's method
restime <- cor.mtest(itta, conf.level = .95, na.rm = TRUE)

#Plot with significance showing - x over NOT significant
png(filename = paste0("Sig_allANT_Correlogram.tif"), width = 10, height = 9.5, units = "in", res = 1200)

print(corrplot.mixed(cor(itta, use = "complete.obs"), p.mat = as.matrix(restime$p), sig.level = .05, tl.col="black", order = "FPC", 
    tl.pos = "lt", tl.cex = 0.5, number.cex = 0.7))

graphics.off()

#Plot without significance showing - only a few not significant because of high n. the significance test is basically pointless. 
#note, the above plot is in a different order than this one because of the significance not working with alphabetical order
png(filename = paste0("allANT_Correlogram.tif"), width = 10, height = 9.5, units = "in", res = 1200)

print(corrplot.mixed(cor(itta, use = "complete.obs"), tl.col="black", order = "alphabet", 
                     tl.pos = "lt", tl.cex = 0.7, number.cex = 0.75))

graphics.off()

##### isolate seropositive data!! :) 

###Seropositivity Thresholds!!!###

#DO NOT exclude any samples or control samples until AFTER seropositivity calculations!! 
#Then can subset antigens, samples, etc from the seropositivity matrix and the final data frame.

#Do seropositivity calculations on NMallsampdata, which has excluded samples removed, but still includes controls.

#need to remove the antigens we aren't using again, from all samp data
bitta <- as.data.frame(tNMallsampdata[,!colnames(tNMallsampdata) %in% rmant])

#confirm that all antigens are in the same order - all TRUE
rownames(finalcut) == colnames(bitta)

#Apply the final cutoff to ALL antigens, ALL dilutions, for test and control samples
SP_all.df <- t(apply(bitta, 1, function(x) ((x > c(finalcut))+0)))

#Make a new data frame where seropositive values will be the data and otherwise it will be NA
SP_all_data.df <- data.frame(matrix(NA, nrow = nrow(bitta), ncol = ncol(bitta)))
rownames(SP_all_data.df) <- rownames(bitta)
colnames(SP_all_data.df) <- colnames(bitta)

for(b in 1:ncol(bitta)){
  for(a in 1:nrow(bitta)){
    if(SP_all.df[[a,b]]==1){
      SP_all_data.df[[a,b]] <- bitta[[a,b]] 
    }
  }
}
remove(a,b)

########################################
############### Plots!!! ###############
########################################

###1. Plot number of seropositive people for each NM antigen, not including negative controls

#Isolate NM antigens and test samples only for seropositivity matrix
SP_NM_test <- filter(target_SP.df, Category == "non_malarial")
SP_NM_test <- tibble::column_to_rownames(SP_NM_test, var = "Name")
SP_NM_test <- SP_NM_test[,sapply(SP_NM_test, is.numeric)]

#data frame of sums of seropositives for each antigen, sorted highest to lowest
SPpeople <- as.data.frame(as.matrix((sort(rowSums(SP_NM_test), decreasing = TRUE))))
SPpeople <- cbind(Target = rownames(SPpeople), SPpeople)
SPpeople$Target <- as.factor(SPpeople$Target)
#explicitly set factor levels to the correct order
SPpeople$Target <- factor(SPpeople$Target, levels = SPpeople$Target[order(-SPpeople$V1)])

png(filename = paste0(study, "_NM_SPpeople.tif"), width = 8, height = 4.5, units = "in", res = 1200)

ggplot(SPpeople, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + 
  ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) +
  theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black"))

graphics.off()

###2. Boxplot of seropositive data from highest to lowest median value
#Add negative controls plotted separately as a point on the same plot? or a small red line?

#isolate data for NM only (it's already test samples only)
#This data frame only has 1182 out of 1325 test samples now, 
#because the others were not seropositive to any antigens and were therefore NA, which is not numeric
SP_NM_test2 <- filter(target_SP_data.df, Category == "non_malarial")
SP_NM_test2 <- tibble::column_to_rownames(SP_NM_test2, var = "Name")
SP_NM_test2 <- SP_NM_test2[,sapply(SP_NM_test2, is.numeric)]

#melt the SP only data with Na.rm = TRUE to organize for ggplot2
SPdatamelt <- melt(as.matrix(SP_NM_test2), na.rm = TRUE)
colnames(SPdatamelt) <- c("Target", "Sample", "Normalized")

#violin plot
ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_NM_All_SP_data_violin.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_violin() + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black"))

graphics.off()

#boxplot
png(filename = paste0(study, "_NM_All_SP_data_box.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) + ylim(0,8)

graphics.off()

#add negative controls as points to the boxplot?
NM_con_data <- filter(target_SP_con_data.df, Category == "non_malarial")
NM_con_data <- tibble::column_to_rownames(NM_con_data, var = "Name")
neg_samples <-c(grep("Neg", colnames(NM_con_data)))
NM_neg_data <- NM_con_data[,neg_samples]
#There are only two values to add to the plot...one for TT and one for RubIV, both the same neg pool
#Not sure if it's worth the trouble of adding these to the plot or just mentioning in the figure caption and text.

png(filename = paste0(study, "_NM_All_SP_data_box_Neg.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) + ylim(0,8)

graphics.off()


