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

library(Amelia)
library(mlbench)
library(caTools)
library(ROCR)

source('~/Desktop/IV Project/CrossVal_Lindsey/ROC_crossval_functions_KG.R')

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

#change the names of the antigens so they are better for everything
mal.all <- tNMallsampdata
colnames(mal.all) <- make.names(colnames(mal.all))

#generate a plot to look for missing data quickly - this sorts columns by most missing data
#as expected we don't have missing data because this study was done with reps = 1
missmap(as.data.frame(mal.all), col=c("blue", "red"), legend=FALSE)

#need to merge this data back with the sample metadata 
mal.meta <- merge(sample_meta_f.df, mal.all, by.x = "sample_id_unique", by.y = "row.names")

### Logistic regression - single antigen - use Etramp5 antigen 1 (Etramp 5 Ag 1 -> Etramp.5.Ag.1)
#the outcome that we are attempting to predict is the column "pcr"

levels(as.factor(mal.meta$pcr))
# "falciparum" "negative" -> the two responses  

#outcome variable must be a factor
mal.meta$pcr <- as.factor(mal.meta$pcr)

### Logistic regression - single antigen - use Etramp5 antigen 1 (Etramp 5 Ag 1 -> Etramp.5.Ag.1)
etramp51.fit <- glm(pcr ~ Etramp.5.Ag.1 , data = mal.meta, family = binomial)

two <- glm(pcr ~ Etramp.5.Ag.1 + HSP40.Ag.1, data = mal.meta, family = binomial)

summary(etramp51.fit)
summary(two)
#the logistic regression is significant, 
#but is this because there are very few positives??

### use cross validation!!! 

#adding in cross validation from Lindsey's script she sent me
formula_mfi <- as.formula("pcr ~ Etramp.5.Ag.1 + HSP40.Ag.1")

m_t30_mfi <- model.crossval(mal.meta, mod.fit="logistic", mod.formula=formula_mfi, perc=.3, outcome.name="pcr", plot=F)
#this is causing an error, even after changing to sample_id from subject_id

roc_t30  <- list(m_t30_mfi)
median_roc_t30 <- list()

med <- which(roc_t30[[1]]$auc.vector==median(roc_t30[[1]]$auc.vector[-1]))[1]
median_roc_t30 <- c(roc_t30[[1]]$performance.list[[med]]@x.values,roc_t30[[1]]$performance.list[[med]]@y.values)


## How to quickly do the logistic regression and cross validation 
#for a number of groups of antigens or scenarios?

glm.probs <- predict(etramp51.fit,type = "response")
glm.two <- predict(two,type = "response")

#confusion matrix - don't know if this is working, don't think it is, 
#I think it's just telling me the number of true and false in the data, not the predictions
table(na.omit(mal.meta$pcr), glm.probs > 0.5)

ROCRpred <- prediction(glm.probs, na.omit(mal.meta$pcr))
ROCRperf <- performance(ROCRpred, 'tpr','fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2,1.7))

ROCRpred2 <- prediction(glm.two, na.omit(mal.meta$pcr))
ROCRperf2 <- performance(ROCRpred2, 'tpr','fpr')
plot(ROCRperf2, colorize = TRUE, text.adj = c(-0.2,1.7))


#need to plot multiple ROC curves on one plot.



