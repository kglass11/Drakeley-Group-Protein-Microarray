#Lou Pk HCC Q-dot Study Analysis 
#August 16, 2018

#Run separately for IgG, IgA, and IgM

###Clear the environnment - OR go to Session > Clear Workspace
rm(list=ls())

#install relevant packages if you haven't used them before
# install.packages("corrgram")
# install.packages("corrplot")
# install.packages("ggbeeswarm")

###Load packages needed for this script
require("gtools")

library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(reshape2)
library(dplyr)
library(corrgram)
library(corrplot)
library(ggbeeswarm)

#set working directory
setwd("I:/Drakeley Group/PROTEIN MICROARRAYS/Experiments/230418 Human Pk case-control Qdot/IgM")
getwd()

#Import data from IgG, IgA, or IgM - this script depends on importing many objects from the end of the processing scripts
load(file = "Pk_HCC_analysis_IgM.RData")

#Add time point to the sample metadata
sample_meta_f.df$day <- 0
sample_meta_f.df$day[grep("D7", sample_meta_f.df$sample_id)] <- 7
sample_meta_f.df$day[grep("D28", sample_meta_f.df$sample_id)] <- 28

#make some columns of sample_meta_f character instead of numeric
sample_meta_f.df$day <- as.character(sample_meta_f.df$day)
sample_meta_f.df$slide_no <- as.character(sample_meta_f.df$slide_no)
sample_meta_f.df$block_rep_1 <- as.character(sample_meta_f.df$block_rep_1)
sample_meta_f.df$block_rep_2 <- as.character(sample_meta_f.df$block_rep_2)

########## Additional Standard Plots ##########

###Plot All standards together in ggplot2 - for ALL THREE ISOTYPES 
Ig <- c("IgA", "Std", "IgM")

for(i in 1:length(Ig)){
  
  type = Ig[i]
  
  norm <- norm.matrix[grep(type, row.names(norm.matrix)),]
  pre <- log.cor.matrix[grep(type, row.names(log.cor.matrix)),]
  
  #normalized
  std1melt <- melt(norm, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_norm_1_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std1melt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
          labs(x = "Sample", y = "Log2(MFI Ratio)", title = "Stds Normalized") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
  #not normalized
  std1premelt <- melt(pre, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_pre_1_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std1premelt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
          labs(x = "Sample", y = "Log2(MFI)", title = "Stds Pre-Normalization") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
}

###### Calculating seropositivity ######

#Seropositivity above a threshold of mean of sample specific buffer spots + 3SD. 
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

#Then can apply the norm_sample_cutoff all antigens
seropos.matrix <- t(apply(norm_sub5.df, 1, function(x) ((x > sub_cutoff)+0)))

#Select Pk PCR+ samples from SEROPOSITIVITY MATRIX
seroposT <- t(seropos.matrix)

SP.meta <- merge(sample_meta_f.df, seroposT, by.y = "row.names", by.x = "sample_id_unique", sort = FALSE)
SP.Pk <- filter(SP.meta, pcr == "Pk")
rownames(SP.Pk) <- SP.Pk$sample_id_unique

SP.Pk.data <- SP.Pk[,(ncol(sample_meta_f.df)+1):ncol(SP.Pk)]

#Select Lou antigens 
SP.Pk.meta <- merge(target_meta.df, t(SP.Pk.data), by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
SP.Pk.Lou.meta <- filter(SP.Pk.meta, Source == "Lou")
SP.Pk.Lou <- tibble::column_to_rownames(SP.Pk.Lou.meta, var="Name")
SP.Pk.Lou <- SP.Pk.Lou[,ncol(target_meta.df):ncol(SP.Pk.Lou)]

#Make a data frame with only seropositive data, NA for everything else
#This is for ALL SAMPLES and ANTIGENS in case you want the other data later
onlySP.df <- as.data.frame(matrix(NA, nrow = nrow(norm_sub5.df), ncol = ncol(norm_sub5.df)))
rownames(onlySP.df) <- rownames(norm_sub5.df)
colnames(onlySP.df) <- colnames(norm_sub5.df)

for(i in 1:nrow(norm_sub5.df)){
  for(k in 1:ncol(norm_sub5.df)){
    if (seropos.matrix[i,k] == 1){
      onlySP.df[i,k] <- norm_sub5.df[i,k]
    }
  }
}

#seropositive data only, filtering out Lou antigens
onlySPmeta.df <- merge(target_meta.df, onlySP.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

onlySPmetaLou.df <- filter(onlySPmeta.df, Source == "Lou")
tacos.df <- tibble::column_to_rownames(onlySPmetaLou.df, var="Name")
tacos.df <- tacos.df[,ncol(target_meta.df):ncol(tacos.df)]

#now select only Pk samples
tacos.meta <- merge(sample_meta_f.df, t(tacos.df), by.y = "row.names", by.x = "sample_id_unique", sort = FALSE)
tacos.Pk <- filter(tacos.meta, pcr == "Pk")
rownames(tacos.Pk) <- tacos.Pk$sample_id_unique

tacos.Pk.data <- tacos.Pk[,(ncol(sample_meta_f.df)+1):ncol(tacos.Pk)]

############ Calculations for total number of people reactive to each Pk antigen ##########

####calculate number of seropositive people at each time point (Pk samples only, Lou antigens only)

#total for all time points --> this is the correct factor order for the plot.   
SPpeople <- as.matrix(sort(rowSums(SP.Pk.Lou), decreasing = TRUE))
  SPpeople <- as.data.frame(SPpeople)
  SPpeople <- tibble::rownames_to_column(SPpeople)
  colnames(SPpeople) <- c("Name", "Total SP")
  
#Add columns for the totals for each day for the antigens in that order
  day <- c("D0","D7","D28")
  
  for(i in 1:3){
  d <- day[i]
  SPday <- SP.Pk.Lou[,grep(d, colnames(SP.Pk.Lou))]
  SPdaysum <- as.matrix(rowSums(SPday))
  colnames(SPdaysum) <- d
  SPpeople <- merge(SPpeople, SPdaysum, by.x = "Name", by.y = "row.names", sort = FALSE)
  }
  remove(i)
  
  #Export this matrix for plotting in PRISM
  write.csv(SPpeople, paste0(study, "_SPpeople4histogram.csv"))
  
############ Calculations and plot for antigen breadth #############

  AntB <- matrix(ncol = 2)
  colnames(AntB) <- c("V1", "day")
  
  #calculate sums for sample for each day and bind together
  for(i in 1:3){
    d <- day[i]
    SPday <- SP.Pk.Lou[,grep(d, colnames(SP.Pk.Lou))]
    AntBday <- as.data.frame(as.matrix(colSums(SPday)))
    AntBday$day <- d
    AntB <- rbind(AntB, AntBday)
  }
  remove(i)
  
  AntB <- AntB[2:nrow(AntB),]
  
  #Plot by day
  
  #set factor order for day
  AntB$day <- factor(AntB$day, levels = as.character(c("D0", "D7", "D28")))
  
  png(filename = paste0(study, "_antigen_breadth_V2.tif"), width = 3, height = 3, units = "in", res = 1200)
  
  ggplot(AntB, aes(x=day, y=V1, color = day)) + geom_violin(color = "black") + 
    geom_beeswarm(cex = 3.5) + 
    theme_bw() + labs(x = "Day", y = "Antigen Breadth", title = iso) +  
    theme(text = element_text(size=11)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
    scale_y_continuous(breaks=c(0,2,4,6,8,10)) +
    theme(legend.position="none") +
    scale_x_discrete(labels=c("0", "7", "28"))
    
  graphics.off()
  
########### Antigen Specific, Correlation between Antibody Response and Age #########

#using seropositive data only, plot data for each time point as a different color
  tacos.Pk.data
  
#need to melt the data for ggplot2
  
  


###########
  
  #parasite count is only taken from day 0, but look at it vs reactivity at all time points. 

########### Correlation Between Antigens, Plot and Stats ############
  
#set up data frame for correlation between antigens for the same points
  #isolate data from lou antigens only and Pk PCR+ samples only!!  
  
  #Make another target.df merged data frame for further use with tag-subtracted values and test samples only
  target2.df <- merge(target_meta.df, norm_sub5.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
  
  #Isolate Lou antigens
  Lou.df <- filter(target2.df, Source == "Lou")
  subLou.df <- tibble::column_to_rownames(Lou.df, var="Name")
  subLou.df <- subLou.df[,sapply(subLou.df, is.numeric)]
  
  #transpose Lou antigens
  subLouT.df <- as.data.frame(t(subLou.df))
  subLouT2.df <- tibble::rownames_to_column(subLouT.df, var = "sample_id_unique")
  Louness.df <- merge(sample_meta_f.df,subLouT2.df, sort = FALSE)
  
  PkLouness <- filter(Louness.df, pcr == "Pk")
  rownames(PkLouness) <- PkLouness$sample_id_unique
  
  PkLouT <- PkLouness[,(ncol(sample_meta_f.df)+1):ncol(PkLouness)]
  
#plot with corrplot package - ALL TIME POINTS COMBINED!
  
  #order FPC means first principle component. I found this order looked the best.I do not know how it was calculated.
  png(filename = paste0(study, "_Pk_antigen_Correlogram.tif"), width = 7, height = 6.5, units = "in", res = 1200)
  
  corrplot.mixed(cor(PkLouT, use = "complete.obs"), tl.col="black", order = "FPC", tl.pos = "lt")
  
  graphics.off()
  
#Combining plot with significance test for correlation coefficients
  #alternative = greater means alternative hypothesis is positive association
  #Right now we are using the default, which is either positive or negative association.
  #default method is pearson, may want to use Kendall or Spearman for our data because better for nonparametric data
  #I tried pearson and spearman and I liked pearson better because we have continuous data 
  #and the spearman (which is rank based) made a lot more things significant.
  res1 <- cor.mtest(PkLouT, conf.level = .95, na.rm = TRUE)
  res2 <- cor.mtest(PkLouT, conf.level = .99, na.rm = TRUE)
  
  #Plots with an X over the correlations which are NOT significant
  png(filename = paste0(study, "_Pk_Correlogram_Sig1.tif"), width = 7, height = 6.5, units = "in", res = 1200)
  
  corrplot.mixed(cor(PkLouT, use = "complete.obs"), tl.col="black", order = "FPC", tl.pos = "lt",
      p.mat = res1$p, sig.level = .05)
  
  graphics.off()
  
  #Plots with a star over the ones that are significant. this looks terrible and I don't know how to make it 
  #only do the star on one half
  png(filename = paste0(study, "_Pk_Correlogram_Sig2.tif"), width = 7, height = 6.5, units = "in", res = 1200)
  
  corrplot.mixed(cor(PkLouT, use = "complete.obs"), tl.col="black", order = "FPC", tl.pos = "lt",
                 p.mat = res1$p, insig = "label_sig", pch.col = "black")
  
  graphics.off()
  
  #Plots with a blank for correlations which are not significant
  png(filename = paste0(study, "_Pk_Correlogram_Sig3.tif"), width = 7, height = 6.5, units = "in", res = 1200)
  
  corrplot.mixed(cor(PkLouT, use = "complete.obs"), tl.col="black", order = "FPC", tl.pos = "lt",
                 p.mat = res1$p, sig.level = .05, insig = "")
  
  graphics.off()
  
###Plot the correlogram separated 1 for each time point!!! 
  cday <- c(0,7,28)
  
  for(i in 1:length(cday)){
    
    #isolate relevant data for each time point
    k = cday[i]
    PkTimeK <- filter(PkLouness, day == k)
    rownames(PkTimeK) <- PkTimeK$sample_id_unique
    PkTimeK <- PkTimeK[,(ncol(sample_meta_f.df)+1):ncol(PkTimeK)]
    
    #run correlation significance test with pearson's method
    restime <- cor.mtest(PkTimeK, conf.level = .95, na.rm = TRUE)
    
    #Plot with significance showing - x over NOT significant
    png(filename = paste0(study,"_", k, "Sig_Pk_Correlogram.tif"), width = 7, height = 6.5, units = "in", res = 1200)
    
    print(corrplot.mixed(cor(PkTimeK, use = "complete.obs"), tl.col="black", order = "FPC", tl.pos = "lt",
                   p.mat = restime$p, sig.level = .05))
    
    graphics.off()
    
    #Plot without significance showing - can add stars when making figure in powerpoint?!?
    png(filename = paste0(study,"_", k, "_Pk_Correlogram.tif"), width = 7, height = 6.5, units = "in", res = 1200)
    
    print(corrplot.mixed(cor(PkTimeK, use = "complete.obs"), tl.col="black", order = "FPC", tl.pos = "lt"))
    
    graphics.off()
    
  }
  remove(i)
  
  