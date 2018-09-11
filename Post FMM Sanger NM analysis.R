#Sanger Data Analysis Script
#Katie Glass
#updated: 9/11/18

#####################################
############DATA ANALYSIS############
#####################################

rm(list=ls())

#I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Data Processed
#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens"
#"/Users/Katie/Desktop/R files from work/100817 Sanger/Sanger NM V2"
setwd("I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens/Sanger NM V2")
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

###### Isolate data for non-malarial antigens and test samples only (no controls). 
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

#Merge with target metadata to filter based on expression tag etc.
target.df <- merge(target_meta.df, Ftestdata, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

#isolate data for non-malarial antigens
NMdatameta <- filter(target.df, Category == "non_malarial")
NMdata <- NMdatameta[,(ncol(target_meta.df)+1):ncol(NMdatameta)]
rownames(NMdata) <- NMdatameta$Name

min(NMdata) #-7.300968
max(NMdata) #8.524411

#Transpose data so antigens are columns - it's now a matrix
tNMdata <- t(NMdata)

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
