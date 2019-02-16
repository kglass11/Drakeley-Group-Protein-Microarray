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

#install.packages("ggpubr")

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

  #also do a similar process for controls so that can get info about seropositivity of control samples 
  contdata <- filter(alldata, sample_type == "control")
  subcont <- contdata[,(ncol(sample_meta_f.df) +1):ncol(contdata)]
  rownames(subcont) <- contdata$sample_id

  tsubcont <- as.data.frame(t(subcont))
  cont.target <- merge(target_meta2.df, tsubcont, by.x = "target_id_unique", by.y ="row.names")
  
  fajitas <- filter(cont.target, Category == "non_malarial" | Category == "malarial")

  fajitas1 <- fajitas[,(ncol(target_meta2.df)+1):ncol(fajitas)]
  row.names(fajitas1) <- make.names(fajitas$Name)
  
  
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
    negcutoffs <- read.csv("Keneba_IgG_v3_negcutoffs.csv")
    negcutoffs <- negcutoffs[,2:ncol(negcutoffs)]
    
    poscutoffs <- read.csv("Keneba_IgG_v3_poscutoffs.csv")
    poscutoffs <- poscutoffs[,2:ncol(poscutoffs)]
    
    multcutoffs <- read.csv("Keneba_IgG_v3_multcutoffs.csv")
    multcutoffs <- multcutoffs[,2:ncol(multcutoffs)]
    
    #identify duplicates and decide which strategy is better for that antigen
    #decisions and explanations are in my written notebook.
    negtwo <- c(as.character(negcutoffs$Name), as.character(multcutoffs$Name))
    negdup <- negtwo[duplicated(negtwo)]
    negdup
    
    postwo <- c(as.character(poscutoffs$Name), as.character(multcutoffs$Name))
    posdup <- postwo[duplicated(postwo)]
    posdup
    
    #make a single data frame of antigen name and SP cutoff 
    #for the negative or positive groups, selecting "plus3SD" and "minus3SD"
    #for dual, selecting "cutoff.pos"
    negcutoffs$Name <- as.character(negcutoffs$Name)
    negfinal <- negcutoffs[!(negcutoffs$Name == "MSP2.Dd2"|negcutoffs$Name == "Mtb.Ag85B"),c(1,3)]
    negfinal$method <- "Neg"
    colnames(negfinal) <- c("Name", "cutoff", "method")
    
    poscutoffs$Name <- as.character(poscutoffs$Name)
    posfinal <- poscutoffs[(poscutoffs$Name == "Influenza.A.H3N2"|poscutoffs$Name == "RSV.GG"),c(1,3)]
    posfinal$method <- "Pos"
    colnames(posfinal) <- c("Name", "cutoff", "method")
    
    multcutoffs$Name <- as.character(multcutoffs$Name)
    twofinal <- multcutoffs[!(multcutoffs$Name == "MSP2.CH150.9"),c(1,3)]
    twofinal$method <- "FMM"
    colnames(twofinal) <- c("Name", "cutoff", "method")
    
    #make one data frame of the final cutoffs for each antigen - total = 105, this number matches the antigen key v7
    SPcutfinal <- rbind(negfinal, posfinal, twofinal)
    
    #put in alphabetical order by antigen
    SPcutfinal <- SPcutfinal[order(SPcutfinal$Name),]
    
####### Make and export tables of seroprevalence for each antigen separated by country and by year for The Gambia
    
  #determine seropositivity status for each sample and each antigen 
  
    #order burritos1 for test samples and fajitas1 for control samples
    sortburritos <- burritos1[order(row.names(burritos1)),]
    sortfajitas <- fajitas1[order(row.names(fajitas1)),]
    
    #confirm that all antigens are in the same order - all TRUE
    SPcutfinal$Name == rownames(sortburritos)
    SPcutfinal$Name == rownames(sortfajitas)
    
    #make seropositivity matrix from all data (test and controls)
    #Apply the final cutoff to ALL antigens, ALL dilutions, for test and control samples
    SP_test.df <- as.data.frame(t(apply(sortburritos, 2, function(x) ((x > c(SPcutfinal$cutoff))+0))))
    SP_test2 <- tibble::rownames_to_column(SP_test.df, "sample_id_unique")
      
    SP_cont.df <- as.data.frame(t(apply(sortfajitas, 2, function(x) ((x > c(SPcutfinal$cutoff))+0))))
    SP_cont2 <- tibble::rownames_to_column(SP_cont.df, "sample_id")
    
    #melt SP data frames
    SPtestmelt <- melt(SP_test2, value.name = "seropositive")
    SPcontmelt <- melt(SP_cont2, value.name = "seropositive")
    
    #check number of rows for SPtestmelt same as allburritomelt - TRUE
    nrow(SPtestmelt) == nrow(allburritomelt)
    
    #merge with allburritomelt to get the SP status for each sample and antigen together with everything else
    burritomeltSP <- merge(allburritomelt, SPtestmelt, by = c("variable", "sample_id_unique"), sort = FALSE)
    
    
   
    
    
    
    
 #copied script from sanger post FMM file:   
    ##### isolate seropositive data!! :) 
    
    ###Seropositivity Thresholds!!!###
    
    #DO NOT exclude any samples or control samples until AFTER seropositivity calculations!! 
    #Then can subset antigens, samples, etc from the seropositivity matrix and the final data frame.
    
    #Do seropositivity calculations on NMallsampdata, which has excluded samples removed, but still includes controls.
    
    #need to remove the antigens we aren't using again, from all samp data
    bitta <- as.data.frame(tNMallsampdata[,!colnames(tNMallsampdata) %in% rmant])
    
    #confirm that all antigens are in the same order - all TRUE
    rownames(finalcut) == colnames(bitta)
    
    #Apply the final cutoff to ALL antigens, ALL dilutions, for test and control samples
    SP_all.df <- t(apply(bitta, 1, function(x) ((x > c(finalcut))+0)))
    
    #Make a new data frame where seropositive values will be the data and otherwise it will be NA
    SP_all_data.df <- data.frame(matrix(NA, nrow = nrow(bitta), ncol = ncol(bitta)))
    rownames(SP_all_data.df) <- rownames(bitta)
    colnames(SP_all_data.df) <- colnames(bitta)
    
    for(b in 1:ncol(bitta)){
      for(a in 1:nrow(bitta)){
        if(SP_all.df[[a,b]]==1){
          SP_all_data.df[[a,b]] <- bitta[[a,b]] 
        }
      }
    }
    remove(a,b)
    
    
    
    
    
    
    
    
    
    
    
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
    
    #compare to original data immediately after background correction to see if anything happened during processing
    originalcor <- t(cor.matrix)
    
    origmeta <- merge(sample_meta_f.df, originalcor, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE)
    HBVoriginal <- merge(HBVelisaMartin, origmeta, by = c("sample_id", "year"))
    
    HBVoriginal$Anti.HBs.IU.L. <- as.numeric(as.character(HBVoriginal$Anti.HBs.IU.L.))
    
    png(filename = paste0(study, "_BackCOR_DATA_HBV.sAg_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(HBVoriginal, aes(x = Anti.HBs.IU.L., y = HBVoriginal$"64_HBV-sAg_1")) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      scale_y_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Log2(Background corrected MFI)", title = "HBV.sAg ELISA vs RPPA raw data") +
      #geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",3], color = "red", size=0.5)+
      #geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=10.1, color = "red", size=0.5) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      geom_smooth(method=lm, na.rm=TRUE)
    
    sp + stat_cor(method = "pearson")
    
    graphics.off()
    
    #do the same steps for foreground data only - the most raw data possible
    originalcor2 <- t(fore.matrix)
    
    origmeta2 <- merge(sample_meta_f.df, originalcor2, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE)
    HBVoriginal2 <- merge(HBVelisaMartin, origmeta2, by = c("sample_id", "year"))
    
    HBVoriginal2$Anti.HBs.IU.L. <- as.numeric(as.character(HBVoriginal2$Anti.HBs.IU.L.))
    
    png(filename = paste0(study, "_RAW_DATA_HBV.sAg_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(HBVoriginal2, aes(x = Anti.HBs.IU.L., y = HBVoriginal2$"64_HBV-sAg_1")) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      scale_y_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Log2(Background corrected MFI)", title = "HBV.sAg ELISA vs RPPA raw data") +
      #geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",3], color = "red", size=0.5)+
      #geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=10.1, color = "red", size=0.5) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      geom_smooth(method=lm, na.rm=TRUE)
    
    sp + stat_cor(method = "pearson")
    
    graphics.off()
    
    #find the values for HBV, Normalized Log2(MFI Ratio) for the positive controls 
    hepBstd <- alldata[grep("HepB_std", alldata$sample_id_unique, ignore.case = TRUE),]
    x <- grep("HBV", c(colnames(hepBstd)), ignore.case = TRUE)
    hepBstd[,x] #the value for the hepBstd is 4.422146
    
    hepBpos.sample <- alldata[grep("382M", alldata$sample_id_unique, ignore.case = TRUE),]
    y <- grep("HBV", c(colnames(hepBpos.sample)), ignore.case = TRUE)
    hepBpos.sample[,y] #the values for the positive sample 382M are 3.269630 4.516825 (2012, 2016 respectively)
    
    
    #We also have data in the main spreadsheet from Martin called "Anti-Hbs Titre" 
    #Only the reactive samples are shown. Let's try and see if this correlates better for some reason.
    #in case the data is not matching or something. - Update - it is matching.
    
    png(filename = paste0(study, "SUBSET_HBV.sAg_RPPAvELISA_regline.tif"), width = 3.8, height = 3.6, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = as.numeric(as.character(anti.HBs.titre)), y = HBV.sAg )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Normalized Log2(MFI Ratio)", title = "HBV.sAg ELISA (subset) vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      geom_smooth(method=lm)
    
    sp + stat_cor(method = "pearson")
    
    graphics.off()
    
####### ELISA data sensitivity and specificity calculations
    
    #set up data frame for sensitivity and specificity calculations
    sens.spec <- as.data.frame(matrix(nrow = 5, ncol = 8))
    colnames(sens.spec) <- c("Antigen", "True_Positives", "True_Negatives", "False_Positives", "False_Negatives", "Sensitivity", "Specificity", "Agreement")
    
    #make a data frame of seropositivity (1) and seronegativity (0) for all antigens
    #What do to about the indeterminate values?? I think just don't count them...
    #set indeterminates to NA
  
    
    
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
