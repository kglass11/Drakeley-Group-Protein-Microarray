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

setwd("/Users/Katie/Desktop/R files from work/120718 Monkey_IVTT")
#"I:/Drakeley Group/PROTEIN MICROARRAYS/Experiments/120718 Monkey_IVTT"

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
sampleinfo$MONKEY <- as.character(sampleinfo$MONKEY)

#make a separate metadata file for green monkeys only
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

# merge this with sample_meta to select by time point
reactive_meta <- merge(sampleinfo, TreactiveSP, by.y = "row.names", by.x = "sample_id", sort = FALSE)

#melt this for ggplot2
antigens <- rownames(reactiveSP.df)
reactivemelt <- melt(reactive_meta, measure.vars = antigens, na.rm = TRUE)

#for each inoculation level separately, for each antigen separately, plot each monkey vs time.  





