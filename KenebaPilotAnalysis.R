#Keneba Pilot Study Data Analysis 
# Nov. 13 2018

#Post QC 

rm(list=ls())

require("gtools")

library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(reshape2)
library(outliers)

workdir <- "/Users/Katie/Desktop/R files from work/Keneba pilot results/IgG_594"
setwd(workdir)

load(file = "KenebaPi_IgG_AfterProcessing.RData")

#Serum Dilution

#Plot all data for each antigen by each dilution

antnames <- rownames(norm_sub4.df)

tnormsub <- as.data.frame(t(norm_sub4.df))

allmeta <- merge(sample_meta_f.df,tnormsub, all.y = TRUE, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE) 






#separate test and controls

