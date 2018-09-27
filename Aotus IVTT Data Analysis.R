#Aotus IVTT Study Data Analysis 

#Katie Glass
#Sept 27, 2018

###Clear the environnment - OR go to Session > Clear Workspace
rm(list=ls())

###Load packages needed for this script
require("gtools")

library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(reshape2)
#library(outliers)
library(plyr)
library(dplyr)
#library(rowr)
library(corrgram)
library(corrplot)
library(ggbeeswarm)

setwd("I:/Drakeley Group/PROTEIN MICROARRAYS/Experiments/120718 Monkey_IVTT")

load("Macaque_IVTT_AfterProcessing.RData")

sample_meta1 <- read.csv("Pvivax_Aotus_Repeated_Expt_Samples_List_092518.csv")

#use this data frame for everything where you need extract metadata
sampleinfo <- merge(sample_meta_f.df, sample_meta1, by.x = "sample_id", by.y = "SAMPLE_ID", all.x = TRUE)

duplicatemeta <- sampleinfo[duplicated(sampleinfo$sample_id),] #there are no duplicated sample IDS

distinct(as.data.frame(sampleinfo$MONKEY)) #there are 12 monkeys, but only numbers for 11

#There are 8 sample IDS for which there is no monkey number or any other metadata 
missingmeta <- filter(sampleinfo, is.na(sampleinfo$MONKEY) & sample_type == "test") 

#make some columns of sampleinfo character instead of numeric
sampleinfo$DAY <- as.character(sampleinfo$DAY)
sampleinfo$slide_no <- as.character(sampleinfo$slide_no)
sampleinfo$block_rep_1 <- as.character(sampleinfo$block_rep_1)




