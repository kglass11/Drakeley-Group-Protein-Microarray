### Keneba Big Study Analysis 
#Katie Glass
#updated: Jan 15 2019

setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis")

#load IgG or IgM 
load("Keneba_IgG_v3_AfterProcessing.RData")

#install latest version of ggpubr from github - this did not work
#gave an error "lazy-load database '/Users/Katie/Library/R/3.5/library/ggpubr/R/ggpubr.rdb' is corrupt"
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

install.packages("ggpubr")

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

library(gridExtra)
library(ggbeeswarm)
library(ggpubr)



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
  avgCV.3 <- mean(CV.3, na.rm = TRUE ) * 100 #12.1552% CV IgG, 11.37715% for IgM

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

###### Prepare the data!!!! Include correct samples etc

  #prepare data so that only test samples are included (Keneba and PHE)
  alldata <- merge(sample_meta_f.df, trans.norm.df, by = "sample_id_unique", all.y = TRUE, sort = FALSE)
  testdata <- filter(alldata, sample_type == "test")
  subtacos <- testdata[,(ncol(sample_meta_f.df) +1):ncol(testdata)]
  rownames(subtacos) <- testdata$sample_id_unique
  
  #Are there any other samples we should remove?!?! what about the 9 duplicates, which one do we use?
  #yes remove these! Only save the first one? or replace with a mean?
  NKC821Ynum <- grep("NKC821Y", rownames(subtacos))
  NKC821Ysamp <- subtacos[NKC821Ynum]
  
  NKC821Ymean <- colMeans(NKC821Ysamp, na.rm = TRUE)
  
  subtacos1 <- subtacos[-NKC821Ynum[2:9],]
  subtacos1[1,] <- NKC821Ymean
  
  #prepare data so that only antigens of interest are included
  subtacosT <- as.data.frame(t(subtacos1))
  targetdata <- merge(target_meta2.df, subtacosT, by.x = "target_id_unique", by.y ="row.names")
  
  burritos <- filter(targetdata, Category == "non_malarial" | Category == "malarial")
  
  #use original antigen names
  #but make antigen names compatible with file export and other stuff in R.
  burritos1 <- burritos[,(ncol(target_meta2.df)+1):ncol(burritos)]
  row.names(burritos1) <- make.names(burritos$Name)
  
  burritosT <- as.data.frame(t(burritos1))
  allburritos <- merge(sample_meta_f.df, burritosT, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE)

######Histograms of all the data for each antigen (includes negative values)  
  
#plot histograms of all data for each antigen separately
  antnames <- colnames(burritosT)
  
  burritomelt <- melt(burritosT)
  allburritomelt <- melt(allburritos, measure.vars = antnames)
  
  setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis/Histograms")
  
  for(i in 1:length(antnames)){
    
    antigen = antnames[i]
    
    ant1 <- filter(allburritomelt, variable == antigen)
    
    plot1 <- ggplot(ant1, aes(x = value)) + geom_histogram(bins = 50, color="black", fill="light blue") + 
                     theme_bw() + 
                     theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
                     labs(x = "Normalized Log2(MFI Ratio)", y = "Count", title = antigen)
    
    plot2 <- ggplot(ant1, aes(x = value, color = Country, fill = Country)) + geom_density(adjust = 3/4, alpha=.25) +
                     theme_bw() + 
                     theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
                     labs(x = "Normalized Log2(MFI Ratio)", y = "Density", title = antigen)
    
    png(filename = paste0(study, "_", antigen,"_histogram_density.tif"), width = 10, height = 3.5, units = "in", res = 1200)
    
    print(grid.arrange(plot1, plot2, ncol = 2, nrow = 1))
    
    graphics.off()
    
    #alternative histogram with percentage displayed on y axis 
    #png(filename = paste0(study, "_", antigen,"_histogram_percent.tif"), width = 3.5, height = 3, units = "in", res = 1200)
    
    #print(ggplot(ant1, aes(x = value)) + geom_histogram(bins = 75, aes(y = (..count..)/sum(..count..)), color="black", fill="light blue") + 
    #labs(x = "Normalized Log2(MFI Ratio)", y = "Percentage of Samples") + scale_y_continuous(labels = scales::percent))
    
    #graphics.off()
    
  }
  
  
####### Correlations of all data between All antigens 
  
  setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis")
  
  #separate into malarial and non-malarial
  mal <- filter(targetdata, Category == "malarial")
  nonmal <- filter(targetdata, Category == "non_malarial")
  
  ###non-malarial correlations
    nonmal1 <- nonmal[,(ncol(target_meta2.df)+1):ncol(nonmal)]
    row.names(nonmal1) <- make.names(nonmal$Name) 
    nonmalT <- as.data.frame(t(nonmal1))
  
    #run correlation significance test with pearson's method
    nonmalcorp <- cor.mtest(nonmalT, conf.level = .95, na.rm = TRUE)
  
    #Plot with significance showing - x over NOT significant
    png(filename = paste0(study, "_Sig_nonmal_Correlogram.tif"), width = 10, height = 9.5, units = "in", res = 1200)
  
    print(corrplot.mixed(cor(nonmalT, use = "complete.obs"), p.mat = as.matrix(nonmalcorp$p), sig.level = .05, tl.col="black", order = "FPC", 
                       tl.pos = "lt", tl.cex = 0.5, number.cex = 0.5))
  
    graphics.off()
  
    #Plot without significance in alphabetical order
    png(filename = paste0(study,"_nonmal_Correlogram.tif"), width = 14, height = 15, units = "in", res = 1200)
  
    print(corrplot.mixed(cor(nonmalT, use = "complete.obs"), tl.col="black", order = "alphabet", 
                       tl.pos = "lt", tl.cex = 0.6, number.cex = 0.4))
  
    graphics.off()
  
  ###malarial correlations
    mal1 <- mal[,(ncol(target_meta2.df)+1):ncol(mal)]
    row.names(mal1) <- make.names(mal$Name) 
    malT <- as.data.frame(t(mal1))
  
    #run correlation significance test with pearson's method
    restime <- cor.mtest(malT, conf.level = .95, na.rm = TRUE)
  
    #Plot with significance showing - x over NOT significant
    png(filename = paste0(study, "_Sig_mal_Correlogram.tif"), width = 10, height = 9.5, units = "in", res = 1200)
  
    print(corrplot.mixed(cor(malT, use = "complete.obs"), p.mat = as.matrix(restime$p), sig.level = .05, tl.col="black", order = "FPC", 
                       tl.pos = "lt", tl.cex = 0.5, number.cex = 0.5))
  
    graphics.off()
  
    #Plot without significance in alphabetical order
    png(filename = paste0(study, "_Mal_Correlogram.tif"), width = 10, height = 9.5, units = "in", res = 1200)
  
    print(corrplot.mixed(cor(malT, use = "complete.obs"), tl.col="black", order = "alphabet", 
                       tl.pos = "lt", tl.cex = 0.6, number.cex = 0.5))
  
    graphics.off()
  
  ###one with luminex Pf antigens only?
    
####### Read in and finalize seropositivity cutoffs for all antigens 
    
    #This information is coming from the R notebook files. 
    #There are 3 different data frames to import, negative, positive, and multiple populations.
    negcutoffs <- load("Keneba_IgG_v3_negcutoffs")
    poscutoffs <- load("Keneba_IgG_v3_poscutoffs")
    #the loading above is not working, need to go back to original files
    #maybe just export .csv files :/
    
    multcutoffs <- read.csv("Keneba_IgG_v3_multcutoffs.csv")
    multcutoffs <- multcutoffs[,2:ncol(multcutoffs)]
    
####### ELISA data - Correlations vs RPPA and Sensitivity and Specificity Calculations
    
    #read in data from .csv files prepared from the excel files Martin sent.
    Martin2012data <- read.csv(file = "2012elisaNKcelldataMartin.csv")
    
    #isolate data from our study which has the matching sample IDs from 2012
    elisasubmelt <- allburritomelt[allburritomelt$year == "2012" & allburritomelt$sample_id %in% Martin2012data$sample_id,]
    
    #Note: there are 191 rows in Martin's data, but only 178 have sample IDs, and only 167
    #of the sample IDs are matching our 2012 sample IDs.... :/
    
    #merge our data with Martin's data by sample_id and year (even though all 2012)
    elisaALLmelt <- merge(elisasubmelt, Martin2012data, by = c("sample_id", "year"))
    
    #actually we wanted our data not melted for the scatter plots
    elisaALL <- merge(allburritos, Martin2012data, by = c("sample_id", "year"))
    
    #Correlation scatter plots using ggplot2, red lines indicate seropositivity, 
    #dashed red line indicates the negative cutoff for RPPA
    #the linear regression line unfortunately removes the 0s, because they are non-finite values on log scale
    #this is an issue for tetanus toxoid
    
#Tetanus Toxoid
    png(filename = paste0(study, "_TT_RPPAvELISA.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = TT.IU.ml, y = TT )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "TT ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=0.01, color = "red", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_TT_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = TT.IU.ml, y = TT )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "TT ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=0.01, color = "red", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      geom_smooth(method=lm)
    
    sp + stat_cor(method = "pearson", label.x = -6, label.y = 8)
    
    graphics.off()
    
#EBV - Nuclear Antigen 1 - 20 negative or 0 values in ELISA data were removed because they cannot be log transformed
    
    elisaALL$EBNA.Titre <- as.numeric(as.character(elisaALL$EBNA.Titre))
    
    png(filename = paste0(study, "_EBV.EBNA.1_RPPAvELISA.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = EBNA.Titre, y = EBV.EBNA.1 )) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "EBV.EBNA.1 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=23, color = "red", size=0.5) +
      geom_vline(xintercept=20, color = "red", linetype = "dashed", size=0.5) +
      ylim(0,10) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_EBV.EBNA.1_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = EBNA.Titre, y = EBV.EBNA.1 )) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA Log2(Titer)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "EBV.EBNA.1 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=23, color = "red", size=0.5) +
      geom_vline(xintercept=20, color = "red", linetype = "dashed", size=0.5) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      ylim(0,10) +
      geom_smooth(method=lm)
      
    sp + stat_cor(method = "pearson")
    
    graphics.off()
    
    #note - when I added ylim(0,8) which puts the axes for the EBNA plot
    #to be the same as the other plots, it is somehow changing the pearson's correlation 
    #and rather than removing 20 rows containing non-finite or missing values
    #it is causing it to remove 54 values, because there are values higher than 8
    #which were being removed. I have set the axis to 10 to fix this problem.
    
    
#CMV - trying CMVpp150 first 
    
    png(filename = paste0(study, "_CMV.pp150_RPPAvELISA.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.pp150 )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "CMV.pp150 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_CMV.pp150_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.pp150 )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "CMV.pp150 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      geom_smooth(method=lm)
    
    sp + stat_cor(method = "pearson")
    
    graphics.off()
    
#CMV again - grade III antigen 
    png(filename = paste0(study, "_CMV.GradeIII_RPPAvELISA.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.GradeIII )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "CMV.GradeIII ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_CMV.GradeIII_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.GradeIII )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "CMV.GradeIII ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      geom_smooth(method=lm)
    
    sp + stat_cor(method = "pearson")
    
    graphics.off()
    
#HBV - surface antigen 
    
    # import the .csv file, this antigen has data from 2012 and 2016 samples
    HBVelisaMartin <- read.csv("HBVelisaMartin.csv")
    
    # merge with allburritos to get data for matching samples only
    elisaHBV <- merge(HBVelisaMartin, allburritos, by = c("sample_id", "year"))
    elisaHBV$Anti.HBs.IU.L. <- as.numeric(as.character(elisaHBV$Anti.HBs.IU.L.))
    
    #plot!
    
    png(filename = paste0(study, "_HBV.sAg_RPPAvELISA.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    ggplot(elisaHBV, aes(x = Anti.HBs.IU.L., y = HBV.sAg)) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "HBV.sAg ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=10.1, color = "red", size=0.5) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_HBV.sAg_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(elisaHBV, aes(x = Anti.HBs.IU.L., y = HBV.sAg)) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "HBV.sAg ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=10.1, color = "red", size=0.5) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      geom_smooth(method=lm, na.rm=TRUE)
      
    sp + stat_cor(method = "pearson", label.x = 3.5, label.y =8)
    
    graphics.off()
    

    
    
    
    #code below not working - source STHDA correlation analysis, correlation of two variables
    #uses package ggpubr
    
    ggscatter(elisaALL, x = EBNA.Titre, y = EBV.EBNA.1, 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson")    
    
    
####### Save the output of the analysis so far
  #save.image(file = "KenebaAnalysis_IgM_v1.RData")
  #save.image(file = "KenebaAnalysis_IgG_v1.RData")
  save.image(file = "KenebaAnalysis_IgG_v2.RData")
    
    
    
    
    
    
    
    
  
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
