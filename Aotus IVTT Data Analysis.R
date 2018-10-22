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

#library(broom)
#library(multcomp)
#library(pbkrtest)
#library(phia)
#library(lmerTest)
#library(multcompView)
#library(MuMIn)

library(afex)
library(caTools)

library(multcompView, PMCMRplus, rcompanion)


# if(!require(psych)){install.packages("psych")}
# if(!require(FSA)){install.packages("FSA")}
# if(!require(lattice)){install.packages("lattice")}
# if(!require(BSDA)){install.packages("BSDA")}
# if(!require(multcompView)){install.packages("multcompView")}
# 
# if(!require(PMCMRplus)){install.packages("PMCMR")}
# if(!require(rcompanion)){install.packages("rcompanion")}

setwd("/Users/Katie/Desktop/R files from work/120718 Monkey_IVTT")
#"I:/Drakeley Group/PROTEIN MICROARRAYS/Experiments/120718 Monkey_IVTT"

load("Macaque_IVTT_AfterProcessing.RData")

#change the study name 
study <- "Aotus_IVTT"

sample_meta1 <- read.csv("Pvivax_Aotus_Repeated_Expt_Samples_List_092518.csv")

#use this data frame for everything where you need extract metadata
sampleinfo1 <- merge(sample_meta_f.df, sample_meta1, by.x = "sample_id", by.y = "SAMPLE_ID", all.x = TRUE)

duplicatemeta <- sampleinfo1[duplicated(sampleinfo1$sample_id),] #there are no duplicated sample IDS

distinct(as.data.frame(sampleinfo1$MONKEY)) #there are 12 monkeys, but only numbers for 11

#There are 8 sample IDS for which there is no monkey number or any other metadata 
missingmeta <- filter(sampleinfo1, is.na(sampleinfo1$MONKEY) & sample_type == "test") 

#it turns out there are actually a lot of duplicate samples, they just have different sample IDs

#isolate duplicate data 
rep1m <- sampleinfo1[duplicated(sampleinfo1[,c("DAY", "INOC_LEVEL", "MONKEY")]),]
rep2m <- sampleinfo1[duplicated(sampleinfo1[,c("DAY", "INOC_LEVEL", "MONKEY")], fromLast = TRUE),]

#it turns out that except for one sample, the duplicates are one serum and one plasma
#rep1m (serum) and rep2m (plasma)

#the 8 samples with no metadata are also being counted as duplicates because all the values are NA

#that one sample is 12824 - for 27050, where it looks like at the day post 275, 
#the inoculation was supposed to be 4,but is written down as 2. I am manually fixing this in the metadata csv file

#Remove duplicate entries automatically (both duplicates) to check if the rest of the samples are serum or plasma
norepinfo <- sampleinfo1[!(duplicated(sampleinfo1[,c("DAY", "INOC_LEVEL", "MONKEY")]) | 
                  duplicated(sampleinfo1[,c("DAY", "INOC_LEVEL", "MONKEY")], fromLast = TRUE)),]

#How many samples are serum and how many are plasma?
length(which(norepinfo$TYPE == "Serum"))
18/(18+70) #0.2045455
length(which(norepinfo$TYPE == "Plasma"))
1 - 18/(18+70) #0.7954545

#bind together plasma duplicate and the unchanged non-duplicates 
#call the final metadata data frame sampleinfo
sampleinfo <- rbind(norepinfo, rep2m)


#This commented section is the script for averaging duplicate spots for APAC, which is similar to this issue
# #identify list of duplicate 1 and duplicate 2, use that to separate the data for rep1 and rep2
# maltargets <- filter(target_meta2.df, Category == "malaria")
# 
# rep1m <- maltargets[duplicated(maltargets$Name),]
# rep2m <- maltargets[duplicated(maltargets$Name, fromLast = TRUE),]
# 
# rep1unique <- rep1m$target_id_unique
# rep2unique <- rep2m$target_id_unique
# 
# rep1data <- norm.matrix[(rownames(norm.matrix) %in% rep1unique),]
# rep2data <- norm.matrix[(rownames(norm.matrix) %in% rep2unique),]
# 
# #get these matrices in the order where can just average the same element on each matrix
# rep1d <- merge(annotation_targets.df, rep1data, by = "row.names",sort = FALSE)
# rep1d <- rep1d[order(rep1d$ID),]
# row.names(rep1d) <- rep1d$target_id_unique
# rep1d <- rep1d[,8:ncol(rep1d)]
# 
# rep2d <- merge(annotation_targets.df, rep2data, by = "row.names",sort = FALSE)
# rep2d <- rep2d[order(rep2d$ID),]
# row.names(rep2d) <- rep2d$target_id_unique
# rep2d <- rep2d[,8:ncol(rep2d)]
# 
# rep1 <- as.matrix(rep1d)
# rep2 <- as.matrix(rep2d)
# 
# ## Calculate correlation coefficient (default is pearson). Deviants are still included.
# repR <- cor(c(rep1), c(rep2), use = "complete.obs")
# print(repR)
# 
# ## Plot replicate 1 v. replicate 2 for each protein or each person and calculate correlation coefficient.
# png(filename = paste0(study, "_replicatescorrelation.tif"), width = 5, height = 4, units = "in", res = 600)
# par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
#     mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
# 
# plot(rep1, rep2, col="red", cex = 0.1)
# mtext(c(paste("Pearson correlation coefficient:", round(repR, digits=4))), side=3, adj=0)
# 
# graphics.off()
# 
# #isolate data that won't be going through the replicate averaging processing:
# norepdata <- norm.matrix[!(rownames(norm.matrix) %in% rep1unique) & !(rownames(norm.matrix) %in% rep2unique),]
# 
# #Average rep1 and rep2
# normaverageI.matrix <- matrix(nrow = nrow(rep1), ncol = ncol(rep1))
# colnames(normaverageI.matrix) = colnames(rep1)
# rownames(normaverageI.matrix) = rownames(rep1)
# 
# normaverageI.matrix <- log2((2^rep1 + 2^rep2)/2)
# #put average data back together with controls/norepdata
# normaverage.matrix <- rbind(normaverageI.matrix, norepdata)
# 
# write.csv(normaverage.matrix, paste0(study, "_average_norm_log_data.csv")) 
# 


#make some columns of sampleinfo character instead of numeric
sampleinfo$DAY <- as.character(sampleinfo$DAY)
sampleinfo$slide_no <- as.character(sampleinfo$slide_no)
sampleinfo$block_rep_1 <- as.character(sampleinfo$block_rep_1)
sampleinfo$MONKEY <- as.character(sampleinfo$MONKEY)

###make a separate metadata file for green monkeys only
mone <- c("30014", "30034", "32028", "32047", "25029", "29012")
mtwo <- c("30014", "30034", "32028", "32047", "25029", "29041")
mthree <- c("30014", "32028", "32047", "25029", "32029")
mfour <- c("30014","27050","32028","32047", "25029","31029")

green1 <- filter(sampleinfo, INOC_LEVEL == 1, 
  (MONKEY == mone[1] | MONKEY == mone[2]| MONKEY == mone[3]| MONKEY == mone[4]| MONKEY == mone[5]| MONKEY == mone[6]))

green2 <- filter(sampleinfo, INOC_LEVEL == 2, 
  (MONKEY == mtwo[1] | MONKEY == mtwo[2]| MONKEY == mtwo[3]| MONKEY == mtwo[4]| MONKEY == mtwo[5]| MONKEY == mtwo[6]))

green3 <- filter(sampleinfo, INOC_LEVEL == 3, 
  (MONKEY == mthree[1] | MONKEY == mthree[2]| MONKEY == mthree[3]| MONKEY == mthree[4]| MONKEY == mthree[5]| MONKEY == mthree[6]))

green4 <- filter(sampleinfo, INOC_LEVEL == 4, 
  (MONKEY == mfour[1] | MONKEY == mfour[2]| MONKEY == mfour[3]| MONKEY == mfour[4]| MONKEY == mfour[5]| MONKEY == mfour[6]))

greenmonkeys <- rbind(green1, green2, green3, green4)

###merge greenmonkey metadata with all antigen data, no cutoffs implemented 

#norm_sub4.df has data with control targets removed, no control samples removed;
#GST subtracted, negative values set to 0 data 
normsub4T <- t(norm_sub4.df)

greenmonkdata <- merge(greenmonkeys, normsub4T, by.y = "row.names", by.x = "sample_id", sort = FALSE) 

###export green monkey data for heatmap in excel 

write.csv(greenmonkdata, file = "greenmonkeydata.csv")

max(greenmonkdata[,16:259]) #5.543199
.5 * max(greenmonkdata[,16:259]) #2.771599

###determine which antigens have at least one value > 0 
isogreen <-  greenmonkdata[,16:259]
rownames(isogreen) <- greenmonkdata$sample_id

greenpos <- apply(isogreen, 1, function(x) ((x > 0)+0))

notzero <- which(rowSums(greenpos) > 0)
length(notzero) #217/244 Pv IVTT antigens reactive in at least 1 sample

nrow(greenmonkdata) * 0.05
#5% of samples = 5.3 samples
#10% of samples = 10.6 samples 

#which antigens are reactive (above 0) in more than 5% of samples (rounding up)? 95 antigens
fiveper <- which(rowSums(greenpos) > 6)
length(fiveper)

#which antigens are reactive (above 0) in more than 10% of samples (rounding up)? 66 antigens
tenper <- which(rowSums(greenpos) > 11)
length(tenper)

#prepare the data to make a heatmap of only the antigens which are reactive in more than 10% of samples
tengreen <- isogreen[,tenper]
mns <- colMeans(tengreen)
tengreen <- tengreen[,order(mns, decreasing = TRUE)]

#prepare data frame to export to excel to make heatmap
subgreenmonk <- merge(greenmonkeys, tengreen, by.y = "row.names", by.x = "sample_id", sort = FALSE) 

###export green monkey data for heatmap in excel 
write.csv(subgreenmonk, file = "subgreenmonk.csv")

###### Plots!!! ######

#for each inoculation level separately, for each antigen separately, plot each monkey vs time.  

#set factor level order for day manually 
subgreenmonk$DAY <- factor(subgreenmonk$DAY, levels = as.character(c("-1", "7", "10", "13", "14", "20", "21", "28", "29", "57")))

#antigens are the top 25 antigens by mean from subgreenmonk / tengreen
antnames <- colnames(tengreen[1:25])

#data is subgreenmonk removing rest of antigens 
subtacos <- subgreenmonk[,1:(16+24)]

#not sure if we need this yet
#subtacos$day <- factor(subtacos$day, levels = as.character(c("0", "7", "28")))

#need to melt the data for ggplot2
subtacosm <- melt(subtacos, measure.vars = antnames, na.rm = TRUE)

#plot antibody response over time, each monkey a different line

for(i in 1:length(antnames)){
  
  antigen = antnames[i]
  
  ant1 <- filter(subtacosm, variable == antigen)
  
  png(filename = paste0(study, "_", antigen,"_time_by_monkey.tif"), width = 8, height = 2.75, units = "in", res = 1200)
  
  print(ggplot(ant1, aes(x = DAY, y = as.numeric(value), color = MONKEY)) +
          geom_point(shape=18, size = 2) +
          geom_line(aes(group = MONKEY)) + 
          facet_wrap(~ INOC_LEVEL, nrow = 1, scales = "free_x") +
          theme_bw() + labs(x = "Day", y = "Log2(MFI Ratio)", title = antigen) + 
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 10, color = "black")) +
          theme(legend.title = element_text(size = 10)) + 
          theme(title = element_text(size = 12, face = "bold")) +
          theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")))
          
  graphics.off()
  
}


######### Seropositivity Cutoffs ########### 

#maybe should look at the data before determining how to do the cutoffs, 
#whether to do them for each sample or each antigen

# seropositivity cutoffs for each antigen based on mean + 3SD of reactivity of uninfected monkeys (time zero)
#use data where the negatives have been set to 0 or not? look at data then decide 

#prepare overall data frame with all data and sample metadata for filtering 
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

#Tailor the norm_sample_cutoff to remove excluded samples
buffer_cutoff.matrix <- as.matrix(norm_sample_cutoff)
rownames(buffer_cutoff.matrix, colnames(norm4.matrix))
sub_cutoff <- buffer_cutoff.matrix[(!rownames(buffer_cutoff.matrix) %in% samples_exclude),]

#Plot the sample cutoffs for samples included in analysis
png(filename = paste0(study, "_Buffer_Cutoffs.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(sub_cutoff, pch='*', col = "blue", ylim=c(0,max(sub_cutoff)*1.25),
     ylab="Seropositivity Cutoff", xlab="Sample (Array)")

graphics.off()

#Then can apply the norm_sample_cutoff all antigens and all samples
seropos.matrix <- t(apply(norm_sub4.df, 1, function(x) ((x > sub_cutoff)+0)))

#At this point, Remove control samples for further analysis
norm_sub6.df <- norm_sub4.df[,colnames(norm_sub4.df) %in% samples_test]
SP.matrix <- seropos.matrix[,1:150]

#based on the plot 3 samples have super low cutoffs, check these out including for antigen breadth
which(sub_cutoff < 0.2) #this is 12 samples
which(sub_cutoff < 0.1) #this is the 3 samples from the plot 12610_125, 12668_141, Neg_1_151
#the only one of the three samples that is included in green monkeys is 12668.
#this sample also had an antigen breadth of zero even though the cutoff is so low. 
#Maybe this sample should be investigated for other QC issues and excluded?
    #buffer (NO DNA) spots -- no issues
    #background - no issues
    #there are no standar plots :( 

###Create a threshold for overall target reactivity
#e.g. To be included in heatmaps and other analyses, perhaps targets should be reacted to by at least 5% of people?
#All Pv antigens
Pv_target_breadth <- rowSums(SP.matrix, na.rm=TRUE)
Pv_target_reactive <- Pv_target_breadth > (ncol(SP.matrix)/100)*5
cat(sum(Pv_target_reactive), "out of", nrow(SP.matrix), "Pv targets are reactive in at least 5% of monkey samples")

#49 out of 244 Pv targets are reactive in at least 5% of monkey samples
#34 out of 244 Pv targets are reactive in at least 10% of monkey samples

#Make a data frame with only seropositive data, NA for everything else
#This is for ALL SAMPLES and ANTIGENS in case you want the other data later
onlySP.df <- as.data.frame(matrix(NA, nrow = nrow(norm_sub4.df), ncol = ncol(norm_sub4.df)))
rownames(onlySP.df) <- rownames(norm_sub4.df)
colnames(onlySP.df) <- colnames(norm_sub4.df)

for(i in 1:nrow(norm_sub4.df)){
  for(k in 1:ncol(norm_sub4.df)){
    if (seropos.matrix[i,k] == 1){
      onlySP.df[i,k] <- norm_sub4.df[i,k]
    }
  }
}

# isolate only SP data for the target reactive only
reactiveSP.df <- onlySP.df[Pv_target_reactive,]
TreactiveSP <- as.data.frame(t(reactiveSP.df))

# merge this with sample_meta for GREEN MONKEYS ONLY to select by time point
reactive_meta <- merge(greenmonkeys, TreactiveSP, by.y = "row.names", by.x = "sample_id", sort = FALSE)

#melt this for ggplot2
antigens <- rownames(reactiveSP.df)
reactivemelt <- melt(reactive_meta, measure.vars = antigens, na.rm = TRUE)

############ Calculations and plot for antigen breadth #############
#taken from Lou HCC Analysis script

#prepare SP.matrix for green monkeys only with metadata
SPt <- as.data.frame(t(SP.matrix))

SP_green <- merge(greenmonkeys, SPt, by.x = "sample_id", by.y = "row.names", sort = FALSE)

#add antigen breadth for each sample as a last column in this data frame

SP_green$Antigen_Breadth <- rowSums(SP_green[(ncol(greenmonkeys)+1):(ncol(greenmonkeys)+ncol(SPt))])

#set factor level order for day manually 
SP_green$DAY <- factor(SP_green$DAY, levels = as.character(c("-1", "7", "10", "13", "14", "20", "21", "28", "29", "57")))

#Plot antigen breadth in a similar style plot to those for ab response over time and inoculation

# labeled by monkey
png(filename = paste0(study, "_Antigen_breadth_monkey.tif"), width = 8, height = 3, units = "in", res = 1200)

print(ggplot(SP_green, aes(x = DAY, y = Antigen_Breadth, color = MONKEY)) +
        geom_point(shape=18, size = 2) +
        geom_line(aes(group = MONKEY)) + 
        facet_wrap(~ INOC_LEVEL, nrow = 1, scales = "free_x") +
        theme_bw() + labs(x = "Day", y = "Number of Reactive Antigens") + 
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 10, color = "black")) +
        theme(legend.title = element_text(size = 10)) + 
        theme(title = element_text(size = 12, face = "bold")) +
        theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")))

graphics.off()

#no different colors 
png(filename = paste0(study, "_Antigen_breadth_Simple.tif"), width = 8, height = 3, units = "in", res = 1200)

print(ggplot(SP_green, aes(x = DAY, y = Antigen_Breadth)) +
        geom_violin(color = "black", scale = "width") + 
        geom_beeswarm(cex = 3.5, color = "red", size = 0.95) +
        facet_wrap(~ INOC_LEVEL, nrow = 1, scales = "free_x") +
        theme_bw() + labs(x = "Day", y = "Number of Reactive Antigens") +  
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
        theme(axis.text = element_text(size = 10, color = "black")) +
        theme(axis.title = element_text(size = 11, face = "bold")) +
        theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")))

graphics.off()

#boxplot plus points
png(filename = paste0(study, "_Antigen_breadth_Box.tif"), width = 8, height = 3, units = "in", res = 1200)

print(ggplot(SP_green, aes(x = DAY, y = Antigen_Breadth)) +
        geom_boxplot(color = "black", outlier.shape = NA) + 
        geom_beeswarm(cex = 3.5, color = "red", size = 0.95) +
        facet_wrap(~ INOC_LEVEL, nrow = 1, scales = "free_x") +
        theme_bw() + labs(x = "Day", y = "Number of Reactive Antigens") +  
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
        theme(axis.text = element_text(size = 10, color = "black")) +
        theme(axis.title = element_text(size = 11, face = "bold")) +
        theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")))

graphics.off()

#boxplot plus points - poster dimensions
png(filename = paste0(study, "_Antigen_breadth_Box_Poster.tif"), width = 9.5, height = 3, units = "in", res = 1200)

print(ggplot(SP_green, aes(x = DAY, y = Antigen_Breadth)) +
        geom_boxplot(color = "black", outlier.shape = NA) + 
        geom_beeswarm(cex = 3.5, color = "red", size = 0.95) +
        facet_wrap(~ INOC_LEVEL, nrow = 1, scales = "free_x") +
        theme_bw() + labs(x = "Day", y = "Number of Reactive Antigens") +  
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
        theme(axis.text = element_text(size = 10, color = "black")) +
        theme(axis.title = element_text(size = 11, face = "bold")) +
        theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")))

graphics.off()

#boxplot no points
png(filename = paste0(study, "_Antigen_breadth_Box_plain.tif"), width = 8, height = 3, units = "in", res = 1200)

print(ggplot(SP_green, aes(x = DAY, y = Antigen_Breadth)) +
        geom_boxplot(color = "black", fill = "light blue") + 
        facet_wrap(~ INOC_LEVEL, nrow = 1, scales = "free_x") +
        theme_bw() + labs(x = "Day", y = "Number of Reactive Antigens") +  
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
        theme(axis.text = element_text(size = 10, color = "black")) +
        theme(axis.title = element_text(size = 11, face = "bold")) +
        theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")))

graphics.off()

####### Statistics #######


# antigen breadth - friedman test with factor of time - separately for each inoculation. 

#filter out the correct data - one inoc level at a time and remove the day = 13 samples

#manual for loop to do each inoculation level - this didn't work as a for loop, got held up 
#on level 2 then level 4

level = 4

Inoc <- filter(SP_green, INOC_LEVEL == level, !DAY == "13")

AntB <- as.data.frame(cbind(Inoc$MONKEY, as.character(Inoc$DAY), Inoc$Antigen_Breadth))
colnames(AntB) <- c("Monkey", "Day", "Antigen_Breadth") 

#make a contingency table to check that there is a unreplicated block design 
AntbX <- xtabs(~ Monkey + Day, data = AntB)
print(AntbX)

AntbX.df <- as.data.frame(AntbX)

#remove monkeys for which there is not a complete design
monkey <- AntbX.df$Monkey[which(AntbX == 0)]

AntB2 <- AntB[!(AntB$Monkey %in% monkey),]

#set factor levels so that the other monkey isn't counted
AntB2$Monkey <- factor(AntB2$Monkey, levels = unique(AntB2$Monkey))

AntbX2 <- xtabs(~ Monkey + Day, data = AntB2)
print(AntbX2)

#run the friedman test
print(friedman.test(Antigen_Breadth ~ Day | Monkey,
              data = AntB2))

#if any friedman tests are significant, can do pairwise comparisons, otherwise move to comparing 
#baseline to each time point only, and use wilcoxan matched pairs test 
#NONE WERE SIGNIFICANT 

#Wilcoxan matched pairs test - day = -1 vs day = other

day = levels(AntB2$Day)

for(i in 2:length(day)){

testday <- day[i]

x = as.numeric(as.character(AntB2$Antigen_Breadth[AntB2$Day == "-1"]))
y = as.numeric(as.character(AntB2$Antigen_Breadth[AntB2$Day == testday]))

print(testday)

print(wilcox.test(x, y, paired = TRUE, alternative = "less"))

}

#end manual for loop

###Stats for the antibody responses for each of the 25 antigens with the highest means. 

##Starting with Friedman test - it turns out that many are significant. Therefore, only 
#doing pairwise comparisons if significant, otherwise not doing anything else for that antigen/inoc combo.

for(k in 1:4){
  
  level = k
  
  Inoc1 <- filter(subtacos, INOC_LEVEL == level, !DAY == "13")
  
  print(level)
  
  for(i in 1:length(antnames)){
    
    antigen = antnames[i]
    print(antigen)
    
    AntD <- as.data.frame(cbind(Inoc1$MONKEY, as.character(Inoc1$DAY), Inoc1[,(15+i)]))
    colnames(AntD) <- c("Monkey", "Day", "Ab_Response") 
    
    AntD$Ab_Response <- as.numeric(as.character(AntD$Ab_Response))
    
    #make a contingency table to check that there is a unreplicated block design 
    AntDX <- xtabs(~ Monkey + Day, data = AntD)
    
    AntDX.df <- as.data.frame(AntDX)
    
    #remove monkeys for which there is not a complete design
    monkey <- AntDX.df$Monkey[which(AntDX == 0)]
    
    AntD2 <- AntD[!(AntD$Monkey %in% monkey),]
    
    #set factor levels so that the other monkey isn't counted
    AntD2$Monkey <- factor(AntD2$Monkey, levels = unique(AntD2$Monkey))
    
    AntDX2 <- xtabs(~ Monkey + Day, data = AntD2)
    
    #run the friedman test
    fried <- friedman.test(Ab_Response ~ Day | Monkey,
                           data = AntD2)
    print(fried)
    
    
    if(is.na(fried$p.value)){} else if (fried$p.value < 0.05){
      PT = pairwiseSignTest(Ab_Response ~ Day, 
                            data   = AntD2,
                            method = "fdr")
      # Adjusts p-values for multiple comparisons;
      # See ?p.adjust for options
      # assumes already ordered by the blocking variable (MONKEY)
      
      print(PT)
      
      #summary of letters - only if any are significant
      if(any(PT$p.adjust <= 0.05)){
      cldList(p.adjust ~ Comparison,
              data = PT,
              threshold  = 0.05)
      }
    }
  }
   
}

remove(i,k, AntD, antigen, AntD2, AntDX, AntDX2, AntDX.df, Inoc1)

### Moving on to wilcoxan matched pairs, only comparing day -1 with the other time points
#similar to antigen breadth, except doing two-sided for everything because we don't know what to expect

for(k in 1:4){
  
  level = k
  
  Inoc1 <- filter(subtacos, INOC_LEVEL == level, !DAY == "13")
  
  print(level)
  
  for(i in 1:length(antnames)){
    
    antigen = antnames[i]
    print(antigen)
    
    AntD <- as.data.frame(cbind(Inoc1$MONKEY, as.character(Inoc1$DAY), Inoc1[,(15+i)]))
    colnames(AntD) <- c("Monkey", "Day", "Ab_Response") 
    
    AntD$Ab_Response <- as.numeric(as.character(AntD$Ab_Response))
    
    #make a contingency table to check that there is a unreplicated block design 
    AntDX <- xtabs(~ Monkey + Day, data = AntD)
    
    AntDX.df <- as.data.frame(AntDX)
    
    #remove monkeys for which there is not a complete design
    monkey <- AntDX.df$Monkey[which(AntDX == 0)]
    
    AntD2 <- AntD[!(AntD$Monkey %in% monkey),]
    
    #set factor levels so that the other monkey isn't counted
    AntD2$Monkey <- factor(AntD2$Monkey, levels = unique(AntD2$Monkey))
    
    AntDX2 <- xtabs(~ Monkey + Day, data = AntD2)
    

    day = levels(AntD2$Day)

    for(i in 2:length(day)){
  
      testday <- day[i]
  
      x = as.numeric(as.character(AntD2$Ab_Response[AntD2$Day == "-1"]))
      y = as.numeric(as.character(AntD2$Ab_Response[AntD2$Day == testday]))
  
      print(testday)
  
      print(wilcox.test(x, y, paired = TRUE, alternative = "two.sided"))
      
    }
  }
  
  
}

### Area Under the Curve of Each Top 25 Antigen - comparing Inoculation 1 and 3.

#only include the 4 monkeys which are the same at inoculation 1 and 3
monks <- c("30014","32028","32047","25029")

#change the order of data in subtacos to be by day
subtacosday <- subtacos[order(subtacos$DAY),]

#isolate inoculation 1 and 3 data for all antigens and time points
#make sure the values for each monkey are in order by day - yes they are finally!
one.auc <- filter(subtacosday, INOC_LEVEL == "1", MONKEY %in% monks)
three.auc <- filter(subtacosday, INOC_LEVEL == "3", MONKEY %in% monks, !DAY == "57")

#make data frames to store the AUC values for each inoculation
AUC1 <- as.data.frame(matrix(nrow=4, ncol = length(antnames)))
colnames(AUC1) <- antnames
rownames(AUC1) <- monks

AUC3 <- as.data.frame(matrix(nrow=4, ncol = length(antnames)))
colnames(AUC3) <- antnames
rownames(AUC3) <- monks

#make a matrix to store the p-values from the paired t-tests
AUCpvals <- matrix(nrow = length(antnames), ncol = 1)
rownames(AUCpvals) <- antnames
colnames(AUCpvals) <- c("p.value")

for(i in 1:length(antnames)){
  antnow <- antnames[i]
  
  #isolate data for each animal and each antigen 
  one.auc1 <- one.auc[one.auc$MONKEY == "30014",colnames(one.auc) %in% antnow]
  one.auc2 <- one.auc[one.auc$MONKEY == "32028",colnames(one.auc) %in% antnow]
  one.auc3 <- one.auc[one.auc$MONKEY == "32047",colnames(one.auc) %in% antnow]
  one.auc4 <- one.auc[one.auc$MONKEY == "25029",colnames(one.auc) %in% antnow]
  
  three.auc1 <- three.auc[three.auc$MONKEY == "30014",colnames(three.auc) %in% antnow]
  three.auc2 <- three.auc[three.auc$MONKEY == "32028",colnames(three.auc) %in% antnow]
  three.auc3 <- three.auc[three.auc$MONKEY == "32047",colnames(three.auc) %in% antnow]
  three.auc4 <- three.auc[three.auc$MONKEY == "25029",colnames(three.auc) %in% antnow]
  
  #use the trapz function, input x and y, get AUC, store data
  AUC1[1,i] <- trapz(x = c(-1,10,14,21,28), y = one.auc1)
  AUC1[2,i] <- trapz(x = c(-1,10,14,21,28), y = one.auc2)
  AUC1[3,i] <- trapz(x = c(-1,10,14,21,28), y = one.auc3)
  AUC1[4,i] <- trapz(x = c(-1,10,13,21,28), y = one.auc4)
  
  AUC3[1,i] <- trapz(x = c(-1,7,14,21,28), y = three.auc1)
  AUC3[2,i] <- trapz(x = c(-1,7,14,21,28), y = three.auc2)
  AUC3[3,i] <- trapz(x = c(-1,7,14,21,28), y = three.auc3)
  AUC3[4,i] <- trapz(x = c(-1,7,14,21,28), y = three.auc4)
  
  #Then use a t test to compare inoculation 1 to 3
  print(antnow)
  AUCtest <- t.test(x = AUC3[,i], y = AUC1[,i], paired = TRUE)
  print(AUCtest)
  
  #store p values from t-tests in the pvals matrix
  AUCpvals[[i]] <- AUCtest$p.value
  
}

#export table of pvalues 
write.csv(AUCpvals, file = "AotusIVTT_AUC_pvalues.csv")

