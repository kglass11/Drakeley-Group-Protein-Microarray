#Sanger Data Analysis - ROC Analysis Test - Malarial Antigens
#Katie Glass
#4/17/19

#####################################
############DATA ANALYSIS############
#####################################

rm(list=ls())

#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens/Sanger NM V2"
#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens"
#
setwd("/Users/Katie/Desktop/R files from work/100817 Sanger/Sanger NM v3")
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
library(ggbeeswarm)

load(file="Sanger.2.Update.RData")

###### Isolate data for malarial antigens and test samples only (no controls) or test samples and controls
###negative values are included here
###this first section is repeated/edited from non-malarial analysis (Post FMM Sanger NM analysis.R)

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

#isolate data for malarial antigens - falciparum only!!! 
#Note, different dilutions are still included.
NMdatameta <- filter(target.df, Plasmodium == "Pf")
NMdata <- NMdatameta[,(ncol(target_meta.df)+1):ncol(NMdatameta)]
rownames(NMdata) <- NMdatameta$Name

min(NMdata) #-8.92014
max(NMdata) #9.999541

#repeat for subdata (test and control samples)
NMallsampmeta <- filter(targetsub.df, Plasmodium == "Pf")
NMallsampdata <- NMallsampmeta[,(ncol(target_meta.df)+1):ncol(NMallsampmeta)]
rownames(NMallsampdata) <- NMallsampmeta$Name

#Transpose data so antigens are columns - it's now a matrix
tNMdata <- t(NMdata)

tNMallsampdata <- t(NMallsampdata)
