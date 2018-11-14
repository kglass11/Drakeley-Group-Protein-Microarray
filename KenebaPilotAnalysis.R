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

workdir <- "/Users/Katie/Desktop/R files from work/Keneba pilot results/Analysis"
setwd(workdir)

#do the analysis for IgG and IgM one at a time. Load either one or the other.
#make sure to include the study name in the filenames of all plots produced.

load(file = "KenebaPi_IgGv2_AfterProcessing.RData")

###### Serum Dilution

#prepare data 
antnames <- rownames(norm_sub4.df)

tnormsub <- as.data.frame(t(norm_sub4.df))

allmeta <- merge(sample_meta_f.df,tnormsub, all.y = TRUE, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE) 

#make sample_dilution a factor
allmeta$sample_dilution <- as.factor(as.character(allmeta$sample_dilution))

#separate test and controls

testmeta <- filter(allmeta, sample_type == "test")

#Plot all data for each antigen by each dilution as a connected line plot

#melt data frame with measure variables as antnames
testmelt <- melt(testmeta, measure.vars = antnames)

#find the max value to use as the ylim of the plots so all are on same axes
max(testmelt$value, na.rm = TRUE) #8.73

#plot all data as a boxplot for each antigen all on one giant plot
png(filename = paste0(study, "_dilution_box.tif"), width = 15, height = 6, units = "in", res = 1200)

print(ggplot(testmelt,  aes(x=variable, y = value, fill = sample_dilution)) + geom_boxplot(outlier.size = 0.2) +
        theme_bw() + labs(x = "Antigen", y = "Log2(MFI Ratio)") + 
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 8, color = "black")) +
        theme(legend.title = element_text(size = 12)) + ylim(0,9) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)))

graphics.off()

#graphs that are antigen specific: 
for(i in 1:length(antnames)){
  
  antigen = antnames[i]
  
  #isolate data for the antigen
  ant1 <- filter(testmelt, variable == antigen)
  
  
}



