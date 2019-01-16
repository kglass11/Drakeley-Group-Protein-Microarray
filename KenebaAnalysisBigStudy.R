### Keneba Big Study Analysis 
#Katie Glass
#updated: Jan 15 2019

setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis")

#load IgG or IgM 
load("Keneba_IgG_v3_AfterProcessing.RData")

#load packages
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(reshape2)
library(outliers)
library(corrgram)
library(corrplot)



#####################################
############DATA ANALYSIS############
#####################################

####### duplicated samples QC
  #there was a replicated sample ID measured 9 times to get the interarray repeatability data 

  #isolate replicated data
  #post normalization and GST subtraction and replicate averaging, nothing set to 0
  sample_rep <- grep("NKC821Y", colnames(quadrulemat))
  NKC821Ydata <- quadrulemat[,sample_rep]

  CV <- apply(NKC821Ydata, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ydata, na.rm = TRUE)
  avgCV <- mean(CV, na.rm = TRUE ) * 100 

  # but it doesn't make sense to get CV of negative numbers or spots where everything is close to 0
  #so let's only use the CV where the mean is > = to 1 (double buffer values)
  #this is getting 60 spots total out of 144
  NKC821Ygr1 <- NKC821Ydata[which(rowMeans(NKC821Ydata, na.rm = TRUE) > 1),]
  CV.3 <- apply(NKC821Ygr1, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ygr1, na.rm = TRUE)
  avgCV.3 <- mean(CV.3, na.rm = TRUE ) * 100 #12.155% CV IgG, 11.38% for IgM

####### IgG vs IgM standard reactivity 

IgG_high_mean <- mean(quadrulemat[grep("IgG Std 1", rownames(quadrulemat)),], na.rm = TRUE)
IgM_high_mean <- mean(quadrulemat[grep("IgM Std 1", rownames(quadrulemat)),], na.rm = TRUE)
IgMfc_high_mean <- mean(quadrulemat[grep("IgM Fc Std 1", rownames(quadrulemat)),], na.rm = TRUE)

if(iso == "IgG"){
  stdratio <- IgG_high_mean/IgM_high_mean
  eratio <- 2^stdratio #7.945517
  stdratioFc <- IgG_high_mean/IgMfc_high_mean
  eratioFc <- 2^stdratioFc #4.89357e+40
}

if(iso == "IgM"){
  stdratio <- IgM_high_mean/IgG_high_mean
  eratio <- 2^stdratio #250270178
  stdratioFc <- IgMfc_high_mean/IgG_high_mean
  eratioFc <- 2^stdratioFc #110338167
}

######Histograms of all the data for each antigen (includes negative values)

  #prepare data so that only test samples are included (Keneba and PHE)
  alldata <- merge(sample_meta_f.df, trans.norm.df, by = "sample_id_unique", all.y = TRUE, sort = FALSE)
  testdata <- filter(alldata, sample_type == "test")
  subtacos <- testdata[,(ncol(sample_meta_f.df) +1):ncol(testdata)]
  rownames(subtacos) <- testdata$sample_id_unique
  
  #Are there any other samples we should remove?!?! what about the 9 duplicates, which one do we use?
  
  #prepare data so that only antigens of interest are included
  subtacosT <- as.data.frame(t(subtacos))
  targetdata <- merge(target_meta2.df, subtacosT, by.x = "target_id_unique", by.y ="row.names")
  
  burritos <- filter(targetdata, Category == "non_malarial" | Category == "malarial")
  
  #use original antigen names
  burritos1 <- burritos[,(ncol(target_meta2.df)+1):ncol(burritos)]
  row.names(burritos1) <- burritos$Name
  
  burritosT <- as.data.frame(t(burritos1))

#plot histograms of all data for each antigen separately
  antnames <- colnames(burritosT)
  
  burritomelt <- melt(burritosT)
  
  for(i in 1:length(antnames)){
    
    antigen = antnames[i]
    
    ant1 <- filter(burritomelt, variable == antigen)
    
    png(filename = paste0(study, "_", antigen,"_histogram.tif"), width = 3.5, height = 3, units = "in", res = 1200)
    
    print(ggplot(ant1, aes(x = value)) + geom_histogram(bins = 75, aes(y = (..count..)/sum(..count..)), color="black", fill="light blue") + 
            labs(x = "Normalized Log2(MFI)", y = "Percentage of Samples") + scale_y_continuous(labels = scales::percent))
          
    graphics.off()
    
    png(filename = paste0(study, "_", antigen,"_histogram2.tif"), width = 3.5, height = 3, units = "in", res = 1200)
    
    print(ggplot(ant1, aes(x = value)) + geom_histogram(bins = 60, color="black", fill="light blue") + 
            labs(x = "Normalized Log2(MFI)", y = "Count"))
    
    graphics.off()
    
  }

#plot density plot with the Gambia and the PHE samples in different colors

#part of Lou macaque script below 
for(i in 1:length(antnames)){
  
  antigen = antnames[i]
  
  ant1 <- filter(subtacosm, variable == antigen)
  
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.PC.tif"), width = 3.5, height = 3, units = "in", res = 1200)
  
  print(ggplot(ant1, aes(x = parasitecount, y = value, color = day)) + geom_point() +
          theme_bw() + labs(x = "Parasite Count at Day 0", y = "Log2(MFI Ratio)", title = antigen) + 
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12))+ xlim(0,7000) + ylim(0,7))
  
  graphics.off()
  
}

#part of sanger script below:
#plot histogram 
png(filename = paste0(study,"_Histogram", antigen, "v1.tif"), width = 4, height = 4, units = "in", res = 600)
par(mfrow=c(1,1))

hist(antibody1,breaks = 40, xlab='Log2(MFI Ratio)',main='')
title(antigen,adj=0.5)
abline(v=fit.ab2$mu, col = "purple", lwd = 1)
abline(v=cutoffSD, col = "blue", lwd = 1.5)
legend("topleft",paste0("cutoff(2SD): ",round(cutoffSD,3)),lty=1,col="blue",cex=0.5,bty="n",y.intersp=1.1,x.intersp=0.2,seg.len=0.5,text.col="blue")

dev.off()

#plot the two clusters density plot
png(filename = paste0(study,"_Kmeans", antigen, "v1.tif"), width = 4, height = 3, units = "in", res = 600)

print(ggplot(clustdata, aes(Ab_Response, ..count.., color = Cluster)) + geom_density(adjust = 3/4) +
        theme_bw() + xlab("Log2(MFI Ratio)") + ggtitle(paste(antigen, round(cut, digits = 3))) +
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        geom_vline(xintercept = cut))

graphics.off()
