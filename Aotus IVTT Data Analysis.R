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

#maybe should look at the data before determining how to do the cutoffs, 
#whether to do them for each sample or each antigen

# seropositivity cutoffs for each antigen based on mean + 3SD of reactivity of uninfected monkeys (time zero)
#use data where the negatives have been set to 0 or not? look at data then decide 

#norm_sub4.df has data with control targets removed, no control samples removed;
#GST subtracted, negative values set to 0 data 

#prepare overall data frame with all data and sample metadata for filtering 
normsub4T <- t(norm_sub4.df)

info_data <- merge(sampleinfo, normsub4T, by.y = "row.names", by.x = "sample_id", sort = FALSE)

#isolate data for time 0 (uninfected monkeys)
uninfected <- filter(info_data, DAY_POST == "-1")
rownames(uninfected) <- uninfected$sample_id_unique

uninf.data <- uninfected[,(ncol(sampleinfo)+1):ncol(uninfected)]

#it won't work to define a cutoff based on this data because there are too many zeros, won't get a good SD. 
#if we had more controls, i.e. controls over all time points, then we could maybe do this.

#going back to the buffer spots method --> in this case will be no DNA control spots, 
#same used for normalizing for Lou's and other studies

#Seropositivity above a threshold of mean of sample specific buffer spots + 3SD. 
#targets_buffer was changed to No DNA spots, so everything with buffer actually refers to that
sample_cutoff <- cor2_buffer_sample_mean + 3*cor2_buffer_sample_sd
log_sample_cutoff <- log2(sample_cutoff)
norm_sample_cutoff <- log_sample_cutoff - log_buffer_sample_mean

#Tailor the norm_sample_cutoff to remove excluded samples and control samples
buffer_cutoff.matrix <- as.matrix(norm_sample_cutoff)
rownames(buffer_cutoff.matrix, colnames(norm4.matrix))
sub_cutoff <- buffer_cutoff.matrix[(!rownames(buffer_cutoff.matrix) %in% samples_exclude),]

#Plot the sample cutoffs for samples included in analysis
png(filename = paste0(study, "_Buffer_Cutoffs.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(sub_cutoff, pch='*', col = "blue", ylim=c(0,max(sub_cutoff)*1.25),
     ylab="Seropositivity Cutoff", xlab="Sample (Array)")

graphics.off()

#### pick up from here tomorrow, the below needs workd

#At this point, Remove control samples for further analysis
norm_sub6.df <- norm_sub5.df[,colnames(norm_sub5.df) %in% samples_test]

#Then can apply the norm_sample_cutoff all antigens
seropos.matrix <- t(apply(norm_sub4.df, 1, function(x) ((x > sub_cutoff)+0)))






