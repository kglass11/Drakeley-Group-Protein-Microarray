### Keneba Big Study Analysis 
#Katie Glass
#updated: Feb 2019

setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis")

#load IgG or IgM 
load("Keneba_IgM_v3_AfterProcessing.RData")

#or load IgG or IgM analysis
#load("KenebaAnalysis_IgG_v3.RData")

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
library(mclust)
library(dbscan)
library(pheatmap)
library(RColorBrewer)



#####################################
############DATA ANALYSIS############
#####################################


#add age bins to sample_meta_f.df - use same age bins as Martin's paper (Blood 2014)
for(i in 1:nrow(sample_meta_f.df)){
  if (is.na(sample_meta_f.df$Age[i])) {sample_meta_f.df$AgeBin[i] <- NA}
  else if (sample_meta_f.df$Age[i] < 3) {sample_meta_f.df$AgeBin[i] <- "1-2"}
  else if (sample_meta_f.df$Age[i] >= 3 & sample_meta_f.df$Age[i] < 6) {sample_meta_f.df$AgeBin[i] <- "3-5"}
  else if (sample_meta_f.df$Age[i] >= 6 & sample_meta_f.df$Age[i] < 10) {sample_meta_f.df$AgeBin[i] <- "6-9"}
  else if (sample_meta_f.df$Age[i] >= 10 & sample_meta_f.df$Age[i] < 13) {sample_meta_f.df$AgeBin[i] <- "10-12"}
  else if (sample_meta_f.df$Age[i] >= 13 & sample_meta_f.df$Age[i] < 16) {sample_meta_f.df$AgeBin[i] <- "13-15"}
  else if (sample_meta_f.df$Age[i] >= 16 & sample_meta_f.df$Age[i] < 20) {sample_meta_f.df$AgeBin[i] <- "16-19"}
  else if (sample_meta_f.df$Age[i] >= 20 & sample_meta_f.df$Age[i] < 26) {sample_meta_f.df$AgeBin[i] <- "20-25"}
  else if (sample_meta_f.df$Age[i] >= 26 & sample_meta_f.df$Age[i] < 40) {sample_meta_f.df$AgeBin[i] <- "26-39"}
  else if (sample_meta_f.df$Age[i] >= 40 & sample_meta_f.df$Age[i] < 55) {sample_meta_f.df$AgeBin[i] <- "40-54"}
  else if (sample_meta_f.df$Age[i] >= 55 & sample_meta_f.df$Age[i] < 75) {sample_meta_f.df$AgeBin[i] <- "55-75"}
}

#explicitly set factor levels for age bins
sample_meta_f.df$AgeBin <- factor(sample_meta_f.df$AgeBin, levels = c("1-2", "3-5","6-9","10-12","13-15","16-19", "20-25","26-39","40-54","55-75"))


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
  
######calculations for sample population table
  
  length(which(sample_meta_f.df$Country=="The Gambia" & sample_meta_f.df$year=="2012" & sample_meta_f.df$Sex == "F" ))
  length(which(sample_meta_f.df$Country=="The Gambia" & sample_meta_f.df$year=="2012" & sample_meta_f.df$Sex == "M" ))
  
  length(which(sample_meta_f.df$Country=="The Gambia" & sample_meta_f.df$year=="2016" & sample_meta_f.df$Sex == "F" ))
  length(which(sample_meta_f.df$Country=="The Gambia" & sample_meta_f.df$year=="2016" & sample_meta_f.df$Sex == "M" ))
  
  length(which(sample_meta_f.df$Country=="The Gambia" & sample_meta_f.df$year=="2012" & sample_meta_f.df$AgeBin == "55-75" ))
  length(which(sample_meta_f.df$Country=="The Gambia" & sample_meta_f.df$year=="2016" & sample_meta_f.df$AgeBin == "55-75"))
  
  
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
    
####### Dengue PCA - dimensions are dengue types 1 - 4 
    #only doing for IgG so far
    #skip for IgM
    
    #isolate dengue data
    denguedata <- burritosT[,grep("DEN", colnames(burritosT))]
    
    #scaling is really important when variables are measured on different scales
    #try with and without scaling for dengue 
    denguePr <- prcomp(na.omit(denguedata), scale = TRUE)
    denguePr
    
    #usually get one fewer components than variables entered, but here got 4 in 4 out
    summary(denguePr)
    
    #repeat without scaling - without scaling is yielding higher proportion of variance in PC1
    #going forward without scaling. PC1 cumulative proportion = 0.9218
    denguePr2 <- prcomp(na.omit(denguedata), scale = FALSE)
    denguePr2
    summary(denguePr2)
    
    #scree plot (sp?) shows variance(square of standard deviation) for each component
    plot(denguePr, type ="l")
    plot(denguePr2, type = "l")
    
    #biplot, which includes eigenvectors for each variable
    biplot(denguePr, scale = 0)
    biplot(denguePr2, scale = 0)
    
    #extract PCA output (item x in the list)
    str(denguePr2)
    
    #only get first two components - 1524 observations were complete (survived na.omit)
    denguedata2 <- cbind(na.omit(denguedata), denguePr2$x[,1:2])
    
    denguemeta <- merge(denguedata2, sample_meta_f.df, by.x = "row.names", by.y = "sample_id_unique", all.x = TRUE)
    
    #this histogram looks really weird and not like how it looked before 
    #:( not sure what's up with this
    hist(denguedata2$PC1)
    
    #make a data frame of just the sample_ids and PC1
    denguePC1 <- as.data.frame(denguedata2$PC1)
    rownames(denguePC1) <- rownames(denguedata2)
    colnames(denguePC1) <- c("DENV.NS1")
    
####### Read in and finalize seropositivity cutoffs for all antigens 

    #This needs to be updated for IgM - not done yet!! 
    if(iso == "IgG"){
    #This information is coming from the R notebook files. 
    #There are 3 different data frames to import, negative, positive, and multiple populations.
    negcutoffs <- read.csv("Keneba_IgG_v3_negcutoffs.csv")
    negcutoffs <- negcutoffs[,2:ncol(negcutoffs)]
    
    poscutoffs <- read.csv("Keneba_IgG_v3_poscutoffs.csv")
    poscutoffs <- poscutoffs[,2:ncol(poscutoffs)]
    
    multcutoffs <- read.csv("Keneba_IgG_v3_multcutoffs.csv")
    multcutoffs <- multcutoffs[,2:ncol(multcutoffs)]
    }
    
    if(iso == "IgM"){
      negcutoffs <- read.csv("Keneba_IgM_v3_negcutoffs.csv")
      negcutoffs <- negcutoffs[,2:ncol(negcutoffs)]
      
      multcutoffs <- read.csv("Keneba_IgM_v3_multcutoffs.csv")
      multcutoffs <- multcutoffs[,2:ncol(multcutoffs)]
    }
    
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

    if(iso == "IgG"){
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
    }
    
    if(iso == "IgM"){
      negcutoffs$Name <- as.character(negcutoffs$Name)
      negfinal <- negcutoffs[,c(1,3)]
      negfinal$method <- "Neg"
      colnames(negfinal) <- c("Name", "cutoff", "method")
      
      multcutoffs$Name <- as.character(multcutoffs$Name)
      twofinal <- multcutoffs[,c(1,3)]
      twofinal$method <- "FMM"
      colnames(twofinal) <- c("Name", "cutoff", "method")
      
      #make one data frame of the final cutoffs for each antigen - total = 105, this number matches the antigen key v7
      SPcutfinal <- rbind(negfinal, twofinal)
    }
    
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
    
    #dup_samples is a data frame from processing that subsetted the samples.df data frame by duplicated samples
    #check whether or not the repeated sample was repeated by year - it was only 2016 so remove it from dup_samples
    dup_samples[dup_samples$sample_id == "NKC821Y",]
    dup_year <- dup_samples[!(dup_samples$sample_id == "NKC821Y"),]
    dup.id <- dup_year$sample_id
    
  #calculate prevalence for each ag as mean of SP data 0s and 1s
    seroprevalence <- as.data.frame(matrix(nrow = nrow(SPcutfinal), ncol = 6))
    colnames(seroprevalence) <- c("Name", "Keneba.2012", "Keneba.2012.paired", "Keneba.2016","Keneba.2016.paired", "England")
  
    #Filter data by country and year
    Keneba.2012 <- filter(burritomeltSP, year == "2012", Country == "The Gambia")
    Keneba.2016 <- filter(burritomeltSP, year == "2016", Country == "The Gambia")
    PHE <- filter(burritomeltSP, Country == "England")
  
    #filter all paired samples
    Keneba.paired <- filter(burritomeltSP, burritomeltSP$sample_id %in% dup.id)
    
    #filter paired by year
    Keneba.2012.paired <- filter(Keneba.paired, year == "2012", Country == "The Gambia")
    Keneba.2016.paired <- filter(Keneba.paired, year == "2016", Country == "The Gambia")
    
    #go through each antigen and calculate the mean of column "seropositive"
    #Keneba.2012
    for(i in 1:nrow(SPcutfinal)){
      antigen = SPcutfinal$Name[i]
      antdata <- filter(Keneba.2012, variable == antigen)
      seroprevalence$Name[i] <- antigen
      seroprevalence$Keneba.2012[i] <- round(mean(antdata$seropositive, na.rm = TRUE)*100, digits = 2)
    }
    
    #Keneba.2016
    for(i in 1:nrow(SPcutfinal)){
      antigen = SPcutfinal$Name[i]
      antdata <- filter(Keneba.2016, variable == antigen)
      seroprevalence$Name[i] <- antigen
      seroprevalence$Keneba.2016[i] <- round(mean(antdata$seropositive, na.rm = TRUE)*100, digits = 2)
    }
    
    #England
    for(i in 1:nrow(SPcutfinal)){
      antigen = SPcutfinal$Name[i]
      antdata <- filter(PHE, variable == antigen)
      seroprevalence$Name[i] <- antigen
      seroprevalence$England[i] <- round(mean(antdata$seropositive, na.rm = TRUE)*100, digits = 2)
    }
    
    #Keneba 2012 paired
    for(i in 1:nrow(SPcutfinal)){
      antigen = SPcutfinal$Name[i]
      antdata <- filter(Keneba.2012.paired, variable == antigen)
      seroprevalence$Name[i] <- antigen
      seroprevalence$Keneba.2012.paired[i] <- round(mean(antdata$seropositive, na.rm = TRUE)*100, digits = 2)
    }
    
    #Keneba 2016 paired 
    for(i in 1:nrow(SPcutfinal)){
      antigen = SPcutfinal$Name[i]
      antdata <- filter(Keneba.2016.paired, variable == antigen)
      seroprevalence$Name[i] <- antigen
      seroprevalence$Keneba.2016.paired[i] <- round(mean(antdata$seropositive, na.rm = TRUE)*100, digits = 2)
    }

    write.csv(seroprevalence, file = paste0(study, "overall.seroprevalence.csv"))
    
    #Get a Keneba overall? How to do this because there are paired inviduals? 
    #Calculate age-adjusted prevalence for Keneba?
    #or are we ok with just the paired samples?
    
    ## Calculate Prevalence by Age Category!!! Add here
    
####### Clustering analysis of ALL data for selected antigens 
    
    #list antigens we are not including in cluster analysis
    #"DENV1.NS1","DENV2.NS1","DENV3.NS1", "DENV4.NS1" --> hold off on removing dengue for now
    rmant <- c("BP.FHA","CMV.Whole","CT.PGP3.D",
               "DT.NIBSC","Etramp.5.Ag.2","Etramp.5.Ag1.GST.var.1","Etramp.5.Ag1.GST.var.2","Etramp.5.Ag1.GST.var.3",
               "Etramp.5.Ag1.His.var.1","Etramp.5.Ag1.His.var.2","Etramp.5.Ag1.His.var.3",
               "EXOB","EXP3","JEV.NS1","MSP11.H103","MSP2.CH150.9","Mtb.TB10.4","Mumps",
               "PHISTC.A","PHISTC.B","Pneumococcal.14","RSV.FG","RSV.GG","SBP1","Tg.GradeIII",
               "Var2CSA","VAR2CSA.DBL5","X800E2A","X800E2C","X800E2D","X800PRISM")
    
    #remove dengue completely
    rmant2 <- c("BP.FHA","CMV.Whole","CT.PGP3.D","DENV1.NS1","DENV2.NS1","DENV3.NS1", "DENV4.NS1",
                         "DT.NIBSC","Etramp.5.Ag.2","Etramp.5.Ag1.GST.var.1","Etramp.5.Ag1.GST.var.2","Etramp.5.Ag1.GST.var.3",
                         "Etramp.5.Ag1.His.var.1","Etramp.5.Ag1.His.var.2","Etramp.5.Ag1.His.var.3",
                         "EXOB","EXP3","JEV.NS1","MSP11.H103","MSP2.CH150.9","Mtb.TB10.4","Mumps",
                         "PHISTC.A","PHISTC.B","Pneumococcal.14","RSV.FG","RSV.GG","SBP1","Tg.GradeIII",
                         "Var2CSA","VAR2CSA.DBL5","X800E2A","X800E2C","X800E2D","X800PRISM")
    
    #prepare the data frame with all of the data - burritosT has the 105 antigens 
    #as columns and all the test samples as rows
    ittacluster <- burritosT[,!colnames(burritosT) %in% rmant]
    nodenguecluster <- burritosT[,!colnames(burritosT) %in% rmant2]
    
    #prepare data for dengue PCA replacement for dengue
    denguePCAcluster <- merge(denguePC1, nodenguecluster, by = "row.names", all.y = TRUE)
    denguePCAcluster <- tibble::column_to_rownames(denguePCAcluster, var = "Row.names")
    
  #### Hierarchical clustering ####
    #need to supply a distance matrix, which is the distance of every point to every other point
    d <- dist(ittacluster)
    
    #usually need to try different algorithms, ward.D2 pre-selected dunno why though
    fitH <- hclust(d, "ward.D2")
    plot(fitH)
    rect.hclust(fitH, k = 3, border = "red")
    
    hclusters <- cutree(fitH, k = 3)
    hclusters
    
  #### model-based clustering ####  
    #this package looks at many models and uses maximizing BIC to select model type and number of clusters
    #requires NA values to be removed...not proceeding with this at this time because
    #the remaining number of complete cases is 414/1540.
    ittamclust <- na.omit(ittacluster)
    
    #fitM <- Mclust(ittamclust)
    #fitM
    #plot(fitM)
    #select plots in the console - see BIC to see how the model was chosen
    
  #### density based clustering ####
     
    #work out input for eps parameter with kNNdist - look for knee / elbow in data
    #look for more directions on choosing these parameters in dbscan info
    #kNNdistplot(ittacluster, k=3)
    #abline(h= 0.7, col = "red", lty = 2)
    #data cannot contain NAs...so can't use this right now! 
    
    #fitD <- dbscan(ittacluster, eps = 0.7 , minPts = 5)
   #fitD
    #noise points here are the points which do not fit in either cluster
    #plot(iris, col = fitD$cluster)
    
    
    #need to see if the clusters mean anything --> do they group with country or age etc
    
    ittaclusterH <- as.data.frame(cbind(hclusters, ittacluster))
    ittaclusterH$hclusters <- as.factor(ittaclusterH$hclusters)
  
    #plot with stat ellipse which shows 95% confidence interval - this is an example of plotting a scatter plot of two antigens
    ggplot(ittaclusterH, aes(x=AMA1, y = CMV.pp150, color = hclusters, fill = hclusters)) + 
      stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
      geom_point(shape = 21, color = "black")
    
    #do other clustering methods, then add all the clusters to ittaclusterH
    
    #then merge with metadata to see if clusters relate to anything
    clustermeta <- merge(sample_meta_f.df, ittaclusterH, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE)
    
    #export clustermeta to make a heatmap in excel
    write.csv(clustermeta, file = paste0(study, "_hclustresults.csv"))
    
    #how do clusters split by country 
    cluster1 <- filter(clustermeta, hclusters == "1")
    cluster2 <- filter(clustermeta, hclusters == "2")
    cluster3 <- filter(clustermeta, hclusters == "3")
    
    #calculate the percentage of each country in each cluster - these results are for without dengue (for IgG)
    #but if I run it again it will be results for with DENV1-4 separately
    length(which(cluster1$Country == "England"))/length(cluster1$Country)*100 #0 IgG, 8.152174 for IgM
    length(which(cluster1$Country == "The Gambia"))/length(cluster1$Country)*100 #100% IgG, 91.84783 IgM
    
    length(which(cluster2$Country == "England"))/length(cluster2$Country)*100 #53.8674, 43.23308 IgM
    length(which(cluster2$Country == "The Gambia"))/length(cluster2$Country)*100 #46.1326, 56.76692 IgM
    
    length(which(cluster3$Country == "England"))/length(cluster3$Country)*100 #2.269044, NA IgM
    length(which(cluster3$Country == "The Gambia"))/length(cluster3$Country)*100 #97.73096, NA IgM
    
  #### Try hclust within pheatmap - make sure to set the method to "ward.D2"
    
    #***** Update ***** Running IgG data with 2 clusters instead of 3 now!! 
    
    #set cutree based on isotype (this is to define the number of clusters)
    if(iso =="IgM"){ncut <- 2}
    if(iso =="IgG"){ncut <- 2}
  
  #first with DENV1-4 included
    pheatmap(ittacluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", scale = "none", cluster_rows = T, 
             cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
             show_colnames = T, na.col = "black", fontsize_col = 8,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmap2.pdf"))
    
    #add annotations on the side for the age category and the country
    #need to get the order of the samples from the pheatmap hclust to get the 
    #matching country and age info
    annotation_colors = list(
      Country = c(England = "lightblue", "The Gambia" = "pink"),
      AgeBin = c("1-2" =  "#9E0142", "3-5" = "#D53E4F","6-9" = "#F46D43",
      "10-12" = "#FDAE61", "13-15" =  "#FEE08B", "16-19" = "#E6F598",
      "20-25" = "#ABDDA4", "26-39" = "#66C2A5","40-54" =  "#3288BD","55-75" ="#5E4FA2"))
         
    heatmapinfo <- pheatmap(ittacluster, silent = TRUE,scale = "none",clustering_method = "ward.D2", 
    clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean",cutree_rows = ncut )
    
    #get clusters from the same hclust within heatmapinfo 
    apple <- as.matrix(cutree(heatmapinfo$tree_row,k=ncut))
    colnames(apple) <- c("cluster")
    
    #this gives a list of lists, can extract the labels from list tree_row to generate annotations data frame
    annotation_labels <- as.data.frame(as.matrix(heatmapinfo$tree_row$labels[heatmapinfo$tree_row$order]))
    colnames(annotation_labels) <- "sample_id_unique"
    
    annotation_info <- merge(annotation_labels, sample_meta_f.df, sort = FALSE, by = "sample_id_unique")
    annotation_info$Country <- as.factor(annotation_info$Country)
    
    #merge with cluster numbers - this will be used for stas on clusters
    cluster_info <- merge(annotation_info, apple, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE)
    cluster_info$cluster <- as.factor(as.character(cluster_info$cluster))
    
    #prepare annotation data frame
    annotation_info_sub <- annotation_info[,c(10,22)]
    rownames(annotation_info_sub) <- annotation_info$sample_id_unique
    
    #save as PDF
    pheatmap(ittacluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", scale = "none", cluster_rows = T, 
             cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
             show_colnames = T, annotation_row = annotation_info_sub, annotation_colors = annotation_colors,
             na.col = "black", fontsize_col = 7 ,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapAnnotated.pdf"))
     
    #save as .tiff for poster - this looks so much better! 
    pheatmap(ittacluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", scale = "none", cluster_rows = T, 
             cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
             show_colnames = T, annotation_row = annotation_info_sub, annotation_colors = annotation_colors,
             na.col = "black", fontsize_col = 7 ,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapAnnotated.tiff"))
    
  #then repeat without dengue at all  
    pheatmap(nodenguecluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", scale = "none", cluster_rows = T, 
             cluster_cols = T, clustering_method = "ward.D2", cutree_rows =ncut,show_rownames = F, 
             show_colnames = T, na.col = "black", fontsize_col = 8,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapNOdenv.pdf"))
    
    #add annotations on the side for the age category and the country
    #need to get the order of the samples from the pheatmap hclust to get the 
    #matching country and age info
    
    heatmapinfo <- pheatmap(nodenguecluster, silent = TRUE,scale = "none",clustering_method = "ward.D2", 
                  clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean",cutree_rows = ncut )
    
    #this gives a list of lists, can extract the labels from list tree_row to generate annotations data frame
    annotation_labels <- as.data.frame(as.matrix(heatmapinfo$tree_row$labels[heatmapinfo$tree_row$order]))
    colnames(annotation_labels) <- "sample_id_unique"
    
    annotation_info <- merge(annotation_labels, sample_meta_f.df, sort = FALSE, by = "sample_id_unique")
    annotation_info$Country <- as.factor(annotation_info$Country)
    
    annotation_info_sub <- annotation_info[,c(10,22)]
    rownames(annotation_info_sub) <- annotation_info$sample_id_unique
    
    pheatmap(nodenguecluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", scale = "none", cluster_rows = T, 
             cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
             show_colnames = T, annotation_row = annotation_info_sub,annotation_colors = annotation_colors,
             na.col = "black", fontsize_col = 7 ,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapNOdenv_Annotated.pdf"))
    
    #Repeat with results from dengue PCA added (denguePCAcluster)  - concluded that this looks terrible
    #the altered scale of dengue PC1 ruins the whole heatmap
      pheatmap(denguePCAcluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", scale = "none", cluster_rows = T, 
             cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
             show_colnames = T, na.col = "black", fontsize_col = 8,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapDENVPC1.pdf"))
    
      #add annotations on the side for the age category and the country
      #need to get the order of the samples from the pheatmap hclust to get the 
      #matching country and age info
    
      heatmapinfo <- pheatmap(denguePCAcluster, silent = TRUE,scale = "none",clustering_method = "ward.D2", 
      clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean",cutree_rows = ncut )
    
      #this gives a list of lists, can extract the labels from list tree_row to generate annotations data frame
      annotation_labels <- as.data.frame(as.matrix(heatmapinfo$tree_row$labels[heatmapinfo$tree_row$order]))
      colnames(annotation_labels) <- "sample_id_unique"
    
      annotation_info <- merge(annotation_labels, sample_meta_f.df, sort = FALSE, by = "sample_id_unique")
      annotation_info$Country <- as.factor(annotation_info$Country)
    
      annotation_info_sub <- annotation_info[,c(10,22)]
      rownames(annotation_info_sub) <- annotation_info$sample_id_unique
    
      pheatmap(denguePCAcluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
             clustering_distance_cols = "euclidean", scale = "none", cluster_rows = T, 
             cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
             show_colnames = T, annotation_row = annotation_info_sub,annotation_colors = annotation_colors,
             na.col = "black", fontsize_col = 7 ,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapDENVPC1_Annotated.pdf"))
    
    #Repeat with SCALING because dengue PCA added (denguePCAcluster) 
      pheatmap(denguePCAcluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
               clustering_distance_cols = "euclidean", scale = "row", cluster_rows = T, 
               cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
               show_colnames = T, na.col = "black", fontsize_col = 8,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapDENVPC1_scale.pdf"))
      
      #add annotations on the side for the age category and the country
      #need to get the order of the samples from the pheatmap hclust to get the 
      #matching country and age info
      
      heatmapinfo <- pheatmap(denguePCAcluster, silent = TRUE,scale = "row",clustering_method = "ward.D2", 
      clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean",cutree_rows = ncut )
      
      #this gives a list of lists, can extract the labels from list tree_row to generate annotations data frame
      annotation_labels <- as.data.frame(as.matrix(heatmapinfo$tree_row$labels[heatmapinfo$tree_row$order]))
      colnames(annotation_labels) <- "sample_id_unique"
      
      annotation_info <- merge(annotation_labels, sample_meta_f.df, sort = FALSE, by = "sample_id_unique")
      annotation_info$Country <- as.factor(annotation_info$Country)
      
      annotation_info_sub <- annotation_info[,c(10,22)]
      rownames(annotation_info_sub) <- annotation_info$sample_id_unique
      
      pheatmap(denguePCAcluster, colors = colors, border_color = NA, clustering_distance_rows = "euclidean", 
               clustering_distance_cols = "euclidean", scale = "row", cluster_rows = T, 
               cluster_cols = T, clustering_method = "ward.D2", cutree_rows = ncut,show_rownames = F, 
               show_colnames = T, annotation_row = annotation_info_sub,annotation_colors = annotation_colors,
              na.col = "black", fontsize_col = 7 ,angle_col = 45,width = 12,filename = paste0("Keneba_",iso, "_pheatmapDENVPC1_scale_Annotated.pdf"))

            
####### Plots of selected epi data vs antibody response - Keneba by year
    
    #Prep Keneba Data
    Kenebamelt <- filter(burritomeltSP, Country == "The Gambia")
    Kenebamelt$year <- as.factor(as.character(Kenebamelt$year))
    
    #Plots for each antigen separately of age vs. antibody response colored by year
    #show SP cutoff on plot 
    
    setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis/IgG Keneba Age v Ab Response")
    
    for(i in 1:length(SPcutfinal$Name)){
      
      antigen = SPcutfinal$Name[i]
      
      #isolate data for the antigen
      ant1 <- filter(Kenebamelt, variable == antigen)
      
      #isolate cutoff for the antigen
      cut1 <- SPcutfinal$cutoff[i]
      
      #age in years (no bins) scatter plot
      png(filename = paste0(study, "_", antigen,"_Ab_vs.age.tif"), width = 3.5, height = 3, units = "in", res = 1200)
      
      print(ggplot(ant1, aes(x = Age , y = value, color = year)) + geom_point(shape = 17, size = 0.75) +
              theme_bw() + labs(x = "Age", y = "Log2(MFI Ratio)", title = paste(antigen, iso)) + 
              theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
              theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
              theme(legend.title = element_text(size = 12))+ 
              #xlim(0,80) + ylim(-2,10) +
              geom_hline(yintercept=cut1, color = "black", size=0.2))
      
      graphics.off()
      
    }
    
    #By gender and year boxplots for all antigens
    
    setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis/IgG Keneba Gender Boxplots")
    
    for(i in 1:length(SPcutfinal$Name)){
      
      antigen = SPcutfinal$Name[i]
      
      #isolate data for the antigen
      ant1 <- filter(Kenebamelt, variable == antigen)
      
      #isolate cutoff for the antigen
      cut1 <- SPcutfinal$cutoff[i]
      
      ant1gender <- filter(ant1, Sex == "M" | Sex == "F")
    
      png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.Gender.tif"), width = 2.7, height = 3, units = "in", res = 1200)
    
      print(ggplot(ant1gender, aes(x = year, y = value, fill = Sex)) + geom_boxplot(outlier.size = 0.4, show.legend=T) +
            theme_bw() + labs(x = "Year", y = "Log2(MFI Ratio)", title = paste(antigen, iso)) + 
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
            theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
            theme(legend.title = element_text(size = 12)) + 
            #ylim(0,8) +
            geom_hline(yintercept=cut1, color = "black", size=0.2))
    
      graphics.off()
    }
    
    #by gender beeswarm and violin plot - this doesn't work by year it looks stupid
    # png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.Gender_V_bee.tif"), width = 3, height = 4, units = "in", res = 1200)
    # 
    # print(ggplot(ant1gender, aes(x = Sex, y = value, color = year)) + geom_violin(scale = "width", color = "black") +
    #         theme_bw() + labs(x = "Year", y = "Log2(MFI Ratio)", title = antigen) + geom_beeswarm(cex = 1, size = 0.5, show.legend = T) +
    #         theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
    #         theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
    #         theme(legend.title = element_text(size = 12)) + ylim(0,8) +
    #         geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
    # 
    # graphics.off()
    
    
    #Antibody response vs Age Category colored by year(2012 or 2016) - the age categories are taken from Martin's paper
    
    setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis/IgG Keneba Age Bin Plots")
    
    for(i in 1:length(SPcutfinal$Name)){
    antigen = SPcutfinal$Name[i]
    
    #isolate data for the antigen
    ant1 <- filter(Kenebamelt, variable == antigen)
    
    #isolate cutoff for the antigen
    cut1 <- SPcutfinal$cutoff[i]

      ant1bin <- filter(ant1, !AgeBin == "NA")
      
      png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.ageBINS_V_bee.tif"), width = 7.5, height = 4, units = "in", res = 1200)
      
      print(ggplot(ant1bin, aes(x = AgeBin, y = value, color = year)) + geom_violin(scale = "width", color = "black") +
              theme_bw() + labs(x = "Age", y = "Log2(MFI Ratio)", title = antigen) + geom_beeswarm(cex = 0.6, size = 0.5) +
              theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
              theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
              theme(legend.title = element_text(size = 12), legend.position="right") + 
              #ylim(0,8) +
              geom_hline(yintercept=cut1, color = "black", size=0.2))
      
      graphics.off()
      
    }
      
    setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis")
      
####### Paired Samples Keneba Data!! :)       
    
    #We will be doing stats on the paired samples data rather than all the data
    #The relevant starting data frames are 
    #Keneba.paired
    #Keneba.2012.paired 
    #Keneba.2016.paired
    
    #choose new age categories for this data and make them refer to 2012 age only 
    #so that the paired samples are always in the same age category
    
    #order both 2012 and 2016 data frames by sample_id and make sure they match
    Keneba.2012.paired <- Keneba.2012.paired[order(Keneba.2012.paired$sample_id),]
    Keneba.2016.paired <- Keneba.2016.paired[order(Keneba.2016.paired$sample_id),]
    
    all(Keneba.2012.paired$sample_id == Keneba.2016.paired$sample_id, na.rm = FALSE)
    #TRUE
    
    #find max age - 49.66
    max(Keneba.2012.paired$Age, na.rm = TRUE)
    
    #new column AgeBin2012
    for(i in 1:nrow(Keneba.2012.paired)){
      if (is.na(Keneba.2012.paired$Age[i])) {Keneba.2012.paired$AgeBin2012[i] <- NA}
      else if (Keneba.2012.paired$Age[i] < 6) {Keneba.2012.paired$AgeBin2012[i] <- "1-5"}
      else if (Keneba.2012.paired$Age[i] >= 6 & Keneba.2012.paired$Age[i] < 13) {Keneba.2012.paired$AgeBin2012[i] <- "6-12"}
      else if (Keneba.2012.paired$Age[i] >= 13 & Keneba.2012.paired$Age[i] < 16) {Keneba.2012.paired$AgeBin2012[i] <- "13-15"}
      else if (Keneba.2012.paired$Age[i] >= 16 & Keneba.2012.paired$Age[i] < 21) {Keneba.2012.paired$AgeBin2012[i] <- "16-20"}
      else if (Keneba.2012.paired$Age[i] >= 21 & Keneba.2012.paired$Age[i] < 36) {Keneba.2012.paired$AgeBin2012[i] <- "21-35"}
      else if (Keneba.2012.paired$Age[i] >= 36 & Keneba.2012.paired$Age[i] < 51) {Keneba.2012.paired$AgeBin2012[i] <- "36-50"}
    }
    
    #explicitly set factor levels for age bins
    Keneba.2012.paired$AgeBin2012 <- factor(Keneba.2012.paired$AgeBin2012, levels = c( "1-5", "6-12","13-15","16-20","21-35","36-50"))
    
    #add same age bins to 2016 data
    Keneba.2016.paired$AgeBin2012 <- Keneba.2012.paired$AgeBin2012
    
    #bind these data frames back together 
    paired.v2 <- rbind(Keneba.2012.paired,Keneba.2016.paired)
    
    #make year a factor instead of numeric
    paired.v2$year <- as.factor(as.character(paired.v2$year))
    
    #Make Boxplots of 2012 vs 2016 paired data by age bin 2012
    setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis/IgG Keneba Paired Data/IgG Age Bin Plots")
    
    for(i in 1:length(SPcutfinal$Name)){
      antigen = SPcutfinal$Name[i]
      
      #isolate data for the antigen
      ant1 <- filter(paired.v2, variable == antigen)
      
      #isolate cutoff for the antigen
      cut1 <- SPcutfinal$cutoff[i]
      
      ant1bin <- filter(ant1, !AgeBin2012 == "NA")
      
      png(filename = paste0(study, "_", antigen,"_Paired_Ab_vs.ageBINS_box.tif"), width = 5.4, height = 3.5, units = "in", res = 1200)
      
      print(ggplot(ant1bin, aes(x = AgeBin2012, y = value, fill = year)) + geom_boxplot(color = "black", outlier.size = 0.4) +
              theme_bw() + labs(x = "Age in 2012", y = "Log2(MFI Ratio)", title = antigen) + 
              #geom_beeswarm(cex = 0.6, size = 0.5) +
              scale_fill_manual(values=alpha(c("blue", "red2"), 0.6)) +
              theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
              theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
              theme(legend.title = element_text(size = 12), legend.position="right") + 
              #ylim(0,8) +
              labs(fill='Year') +
              theme(plot.title = element_text(color = "black", face = "bold"),
              axis.title.x = element_text(color = "black", face = "bold"),
              axis.title.y = element_text(color = "black", face = "bold"))+
              geom_hline(yintercept=cut1, color = "black", size=0.2, linetype = "dashed"))
      
      graphics.off()
      
    }
    
    setwd("/Users/Katie/Desktop/R files from work/Keneba main results/Keneba Analysis")
    

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
    png(filename = paste0(study, "_TT_RPPAvELISA.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = TT.IU.ml, y = TT )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Log2(MFI Ratio)", title = "TT ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=0.01, color = "red", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_TT_RPPAvELISA_regline.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = TT.IU.ml, y = TT )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Log2(MFI Ratio)", title = "TT ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "TT",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=0.01, color = "red", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
      geom_smooth(method=lm)
    
    sp + stat_cor(method = "pearson", label.x = -6, label.y = 8)
    
    graphics.off()
    
#EBV - Nuclear Antigen 1 - 20 negative or 0 values in ELISA data were removed because they cannot be log transformed
    
    elisaALL$EBNA.Titre <- as.numeric(as.character(elisaALL$EBNA.Titre))
    
    png(filename = paste0(study, "_EBV.EBNA.1_RPPAvELISA.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = EBNA.Titre, y = EBV.EBNA.1 )) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(RU/mL)" , y = " RPPA - Log2(MFI Ratio)", title = "EBV.EBNA.1 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=23, color = "red", size=0.5) +
      geom_vline(xintercept=20, color = "red", linetype = "dashed", size=0.5) +
      ylim(0,10) +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) 
 
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_EBV.EBNA.1_RPPAvELISA_regline.tif"),width = 3.23, height = 3.35, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = EBNA.Titre, y = EBV.EBNA.1 )) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA Log2(RU/ML)" , y = " RPPA - Log2(MFI Ratio)", title = "EBV.EBNA.1 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "EBV.EBNA.1",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=23, color = "red", size=0.5) +
      geom_vline(xintercept=20, color = "red", linetype = "dashed", size=0.5) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
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
    
    png(filename = paste0(study, "_CMV.pp150_RPPAvELISA.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.pp150 )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Log2(MFI Ratio)", title = "CMV.pp150 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_CMV.pp150_RPPAvELISA_regline.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.pp150 )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Log2(MFI Ratio)", title = "CMV.pp150 ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.pp150",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
      geom_smooth(method=lm)
    
    sp + stat_cor(method = "pearson")
    
    graphics.off()
    
#CMV again - grade III antigen 
    png(filename = paste0(study, "_CMV.GradeIII_RPPAvELISA.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.GradeIII )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Log2(MFI Ratio)", title = "CMV.GradeIII ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_CMV.GradeIII_RPPAvELISA_regline.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    sp <- ggplot(elisaALL, aes(x = as.numeric(as.character(HCMV.IgG.Titre)), y = CMV.GradeIII )) + geom_point(color = "darkblue", size = 0.7) + 
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(Titer)" , y = " RPPA - Log2(MFI Ratio)", title = "CMV.GradeIII ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "CMV.GradeIII",2], color = "red", linetype = "dashed", size=0.5)+
      #geom_vline(xintercept=?, color = "red", size=0.5) +
      #geom_vline(xintercept=?, color = "red", linetype = "dashed", size=0.5) +
      scale_x_continuous(trans='log2') +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
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
    
    png(filename = paste0(study, "_HBV.sAg_RPPAvELISA.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    ggplot(elisaHBV, aes(x = Anti.HBs.IU.L., y = HBV.sAg)) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Log2(MFI Ratio)", title = "HBV.sAg ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=10.1, color = "red", size=0.5) +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))
    
    graphics.off()
    
    #repeat plot with regression line, 95% confidence intervals, pearson's R, and p value of correlation test
    png(filename = paste0(study, "_HBV.sAg_RPPAvELISA_regline.tif"), width = 3.23, height = 3.35, units = "in", res = 1200)
    
    sp <- ggplot(elisaHBV, aes(x = Anti.HBs.IU.L., y = HBV.sAg)) + geom_point(color = "darkblue", size = 0.7) + 
      scale_x_continuous(trans='log2') +
      theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      labs(x = "ELISA - Log2(IU/mL)" , y = " RPPA - Log2(MFI Ratio)", title = "HBV.sAg ELISA vs RPPA") +
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",3], color = "red", size=0.5)+
      geom_hline(yintercept=multcutoffs[multcutoffs$Name == "HBV.sAg",2], color = "red", linetype = "dashed", size=0.5)+
      geom_vline(xintercept=10.1, color = "red", size=0.5) +
      theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm")) +
      theme(plot.title = element_text(color = "black", face = "bold"),
            axis.title.x = element_text(color = "black", face = "bold"),
            axis.title.y = element_text(color = "black", face = "bold"))+
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
    sens.spec <- as.data.frame(matrix(nrow = 4, ncol = 8))
    colnames(sens.spec) <- c("Antigen", "True_Positives", "True_Negatives", "False_Positives", "False_Negatives", "Sensitivity", "Specificity", "Agreement")
    
    #isolate data to use in the ELISA analysis
    ELISAant <- cbind(elisaALL$CMV.GradeIII, elisaALL$CMV.pp150, elisaALL$EBV.EBNA.1, elisaALL$TT)
    colnames(ELISAant) <- c("CMV.GradeIII","CMV.pp150", "EBV.EBNA.1", "TT")
    
    #isolate SP cutoffs for those antigens
    SPcuts <- multcutoffs[(multcutoffs$Name %in% c("CMV.GradeIII","CMV.pp150", "EBV.EBNA.1", "TT")),3]
    
    #isolate SN cutoffs for those antigens
    SNcuts <- multcutoffs[(multcutoffs$Name %in% c("CMV.GradeIII","CMV.pp150", "EBV.EBNA.1", "TT")),2]
    
    #make SP binary data frame, set indeterminates to NA
    SPbinary <- as.data.frame(matrix(ncol=ncol(ELISAant), nrow = nrow(ELISAant)))
    colnames(SPbinary) <- colnames(ELISAant)
    
    for(i in 1:ncol(ELISAant)){
      for(j in 1:nrow(ELISAant)){
        if(is.na(ELISAant[j,i])){
          SPbinary[j,i] <- NA
        }else if(ELISAant[j,i] >= SPcuts[i]){
          SPbinary[j,i] <- 1
          }else if(ELISAant[j,i] <= SNcuts[i]){
            SPbinary[j,i] <- 0
          }else if(ELISAant[j,i] > SNcuts[i] & ELISAant[j,i] < SPcuts[i]){
            SPbinary[j,i] <- NA
          }
      }
    }
    
    #add the columns for the data from Martin
    SPbinary <- cbind(SPbinary, elisaALL$EBNA.Result,elisaALL$TT.Result)
    SPbinary$CMV.Result <- (!(is.na(elisaALL$HCMV.IgG.Titre)))+0
    
    #fill out sensitivity and specificity data frame
    sens.spec$Antigen <- c("CMV.GradeIII","CMV.pp150", "EBV.EBNA.1", "TT")
    
    #CMV.GradeIII - there are no negatives by ELISA in the data given
    sens.spec$True_Positives[1] <- length((which(SPbinary$CMV.GradeIII == SPbinary$CMV.Result)))
    sens.spec$True_Negatives[1] <- 0
    sens.spec$False_Positives[1] <- 0
    sens.spec$False_Negatives[1] <- length(which(!(SPbinary$CMV.GradeIII == SPbinary$CMV.Result)))
    
    #CMV.pp150 - there are no negatives by ELISA in the data given
    sens.spec$True_Positives[2] <- length((which(SPbinary$CMV.pp150 == SPbinary$CMV.Result)))
    sens.spec$True_Negatives[2] <- 0
    sens.spec$False_Positives[2] <- 0
    sens.spec$False_Negatives[2] <- length(which(!(SPbinary$CMV.pp150 == SPbinary$CMV.Result)))
    
    #EBV.EBNA.1
    sens.spec$True_Positives[3] <- length(which(SPbinary$EBV.EBNA.1 == (SPbinary$`elisaALL$EBNA.Result` == 1)))
    sens.spec$True_Negatives[3] <- length(which(SPbinary$EBV.EBNA.1 == (SPbinary$`elisaALL$EBNA.Result` == 0)))
    sens.spec$False_Positives[3] <- length(which((SPbinary$EBV.EBNA.1 == 1) & (SPbinary$`elisaALL$EBNA.Result` == 0)))
    sens.spec$False_Negatives[3] <- length(which((SPbinary$EBV.EBNA.1 == 0) & (SPbinary$`elisaALL$EBNA.Result` == 1)))
    
    #TT
    sens.spec$True_Positives[4] <- length(which(SPbinary$TT == (SPbinary$`elisaALL$TT.Result` == 1)))
    sens.spec$True_Negatives[4] <- length(which(SPbinary$TT == (SPbinary$`elisaALL$TT.Result` == 0)))
    sens.spec$False_Positives[4] <- length(which((SPbinary$TT == 1) & (SPbinary$`elisaALL$TT.Result` == 0)))
    sens.spec$False_Negatives[4] <- length(which((SPbinary$TT == 0) & (SPbinary$`elisaALL$TT.Result` == 1)))
    
    #calculate the other quantities
    #sensitivity = true positives / (true positives + false negatives)
    #specificity = true negatives / (true negatives + false positives)
    #percent agreement = % (true positives + true negatives) / total sample number (that aren't NA)
    
    for(i in 1:nrow(sens.spec)){
      sens.spec$Sensitivity[i] <- sens.spec$True_Positives[i]/(sens.spec$True_Positives[i] + sens.spec$False_Negatives[i])
      sens.spec$Specificity[i] <- sens.spec$True_Negatives[i]/(sens.spec$True_Negatives[i] + sens.spec$False_Positives[i])
      sens.spec$Agreement[i] <- (sens.spec$True_Positives[i] + sens.spec$True_Negatives[i])/(sum(sens.spec[i,2:5]))*100
    }
      
      
####### Save the output of the analysis so far
  save.image(file = "KenebaAnalysis_IgM_v2.RData")
  #save.image(file = "KenebaAnalysis_IgG_v1.RData")
  #save.image(file = "KenebaAnalysis_IgG_v3.RData")
    
    
  