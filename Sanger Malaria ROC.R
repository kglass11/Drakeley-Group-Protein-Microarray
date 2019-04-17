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