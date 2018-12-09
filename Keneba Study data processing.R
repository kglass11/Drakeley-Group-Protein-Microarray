###Combined script for reading in and processing microarray data to prepare for analysis
  #IgG or IgM
  #Last updated December 4, 2018 for Keneba Big Study

#Updated because the new antibodies have IgG at 594 and IgM at 488

###Create a folder, into which you should copy this script, your .gpr files (named 'Slide 1.gpr', 'Slide 2.gpr' etc.), 
# and your sample list, sample metadata, and target (antigen) metadata csv files.

# Your sample list file needs six columns. Their names must be exactly as written here, though the order of the samples does not matter:
#1.slide_no
#2.sample_id
#3.block_rep_1
#4.block_rep_2
#5.exclude - where samples you want to exclude are labeled "yes"
#6.sample_type - where samples are labeled as "test" or "control"

# In your target metadata file, the targets must be listed in the same order as they are in the .gpr files. 
# If you have two identical blocks printed (for duplicates), only list block 1 in your target metadata file. 

###Clear the environnment - OR go to Session > Clear Workspace
rm(list=ls())

###Install any packages you may need for this script. Go to Tools, install packages, 
#then install from CRAN repository or local file or run:
#install.packages(c("dplyr","gtools","contrast", "beeswarm", "mixtools", "gplots", "ggplot2", "gcookbook", "reshape2"))

#install.packages("outliers")

##Install limma if you haven't already
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")

###Load packages needed for this script
require("gtools")

library(limma)
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

### Define variables based on your study that will be used later in the script
# define working directory character vector, example "I:/Drakeley Group/Protein microarrays/Experiments/030417 Ghanaian samples/RepeatProcessingMay21KG"
workdir <- "/Users/Katie/Desktop/R files from work/Keneba main results"

# define a shorthand name for your study which will be appended in the file name of all exported files
#include isotype in the study name!
study <- "Keneba_IgG.2"

#define isotype
iso <- "IgG"

#define file name for sample IDs character vector, example "Analysis sample list 2.csv"
sample_file <- "Keneba main sample list.csv"

#define file name for sample file + additional metadata (character vector)
meta_file <- "Keneba Sample Metadata.csv"

#define file name for antigen list file with additional info about targets.
target_file <- "KenebaAntigenKeyv4.csv" 

#number of technical replicates for the study (usually 1 or 2)
reps <- 4

#define number of blocks per slide
index_block <- 32


### Set the working directory for the folder where .gpr files are. Can check working
#directory after with getwd()
setwd(workdir)
getwd()

#####################################
###READING IN YOUR MICROARRAY DATA###
#####################################

###Identify all .gpr files in your set path. default path for list.files() is the working directory.
slide_ids <- list.files(pattern="*.gpr")

###Read in the contents of each of the .gpr files you just identified
#This simply combines the data in a list format
#If your gpr files differ in terms of the number of rows of information before the header, you need to edit this command
#i.e skip=34 means it will skip the first 34 rows, and header=T means it will read the 35th ow as a header, and the 36th as the first row of data
# which is what currently happens in our GPR files.
#Note: slide_no_temp should be the last slide added to the list
#Slide number is assigned based on the file name. As files are created without editing by the genepix software - this seems like the easiest way of identifying files.
slides_list <- list()
for(i in 1:length(slide_ids)) { 
  
  slides_list[[i]] <- read.table(slide_ids[i],skip=34,sep="\t",header=T)
  slide_no_temp <- substr(slide_ids[[i]],7,nchar(slide_ids[[i]])-4)
  slides_list[[i]][,"slide_no"] <- slide_no_temp
}
remove(i)

###Name your slides in this list, according to their file names
names(slides_list) <- slide_ids

#This part was for Keneba pilot study: 
#change the column names for slide 4 to be the same as the others because 
#they were messed up in genepix software when re-extracting the data
#colnames(slides_list[[4]]) <- colnames(slides_list[[3]])

###Bind all data from the slide data list (slides.list) into a single dataframe
#you may get a warning after this step, invalid factor level, this is not a problem!
slides_all.df <- c()
for(i in 1:length(slides_list)) { 
  
  slides_all.df <- rbind(slides_all.df,slides_list[[i]])
}
remove(i)

###Read in list of sample IDs
samples1.df <- read.csv(sample_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

###Read in sample metadata file
sample_meta1.df <- read.csv(meta_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

###Read in target metadata file
target_meta1.df <- read.csv(target_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

#remove duplicates from target metadata
target_meta.df <- distinct(target_meta1.df)

###Processing sample list:

##remove extra NAs in the imported sample list 
samples.df <- samples1.df[!is.na(samples1.df$slide_no),]

###Create a vector listing all of your samples, in the order they appear in your samples_list file
#In our samples_list files, the sample ID is always in column two - which is why that column is picked out here
samples <- as.character(samples.df[1:nrow(samples.df), 2])

###Now, make a new sample variable that is unique (to avoid issues with multiple blanks, etc.) 
 # Add this column at the end of the sample list file, not at a prespecified column number
samples.df$sample_id_unique <- c()
samples_unique <- c(paste(as.character(samples.df[1:nrow(samples.df), 2]), rownames(samples.df), sep = "_"))
samples.df$sample_id_unique <- samples_unique

### Cleaning sample metadata (epi data)
#Remove one of the duplicates where the whole row is duplicated (i.e. same thing listed twice)
sample_meta2.df <- distinct(sample_meta1.df)

#do separately for each year because there are duplicate entries for every year
sample_meta_2012 <- filter(sample_meta2.df, year == "2012")
sample_meta_2016 <- filter(sample_meta2.df, year == "2016")
#Export a table of duplicate entries that do not match (i.e. there is a problem with epi data) 
duplicate_metadata <- sample_meta_2012[(duplicated(sample_meta_2012$sample_id)| duplicated(sample_meta_2012$sample_id, fromLast=TRUE)),]
write.csv(duplicate_metadata, file = paste0(study, "_duplicate_metadata2012.csv"))

duplicate_metadata2016 <- sample_meta_2016[(duplicated(sample_meta_2016$sample_id)| duplicated(sample_meta_2016$sample_id, fromLast=TRUE)),]
write.csv(duplicate_metadata2016, file = paste0(study, "_duplicate_metadata2016.csv"))

## there is no duplicate metadata with mismatched info!!! skip the following
#Remove duplicate entries automatically (both duplicates)
#sample_meta3.df <- sample_meta2.df[!(duplicated(sample_meta2.df$sample_id) | duplicated(sample_meta2.df$sample_id, fromLast=TRUE)),]

#Set exclude = "yes" in sample list file for duplicates that don't match in metadata
#dup <- unique(duplicate_metadata$sample_id)
#for(i in 1:length(dup)){
#  samples.df$exclude[which(samples.df$sample_id == dup[i])] <- "yes"
#}
#remove(i)

sample_meta3.df <- sample_meta2.df

### Merge the sample list file and the sample metadata file to include the appropriate metadata
#The duplicate metadata will now be listed as NA, with exclude = yes
sample_meta.df <- merge(samples.df, sample_meta3.df, by = c("sample_id", "year"), all.x = TRUE, sort = FALSE)

###Create vectors indicating the number of slides, blocks, and samples
#Slide and sample number are determined automatically from the data you input, whereas block number is manual in this instance
index_slide <- as.numeric(length(slides_list))
index_sample <- as.numeric(length(samples))

###Assign your sample_ids to each row of the combined slide data (slides_all.df)
#The order of data in your samples.df file is irrelevant, as long as each sample ID is correctly matched to its slide and block numbers
slides_all.df$Sample <- NA

for(i in 1:nrow(slides_all.df)){
  print(i)
  block <- slides_all.df$Block[[i]]
  slide <- as.numeric(slides_all.df$slide_no[[i]])
  temp <- which(samples.df$slide_no == slide & (samples.df$block_rep_1 == block | samples.df$block_rep_2 == block))
  name <- samples.df$sample_id_unique[temp]
  slides_all.df$Sample[[i]] <- name
}
remove(i,block,slide,temp,name)

###Write slides_all.df to a file to keep as a csv in your directory
write.csv(slides_all.df,file=paste0(study,"_slidesall_combinedGPR.csv"), row.names=T)

###Save slides_all.df as an R object to be loaded later so that you don't have to redo that part
save(slides_all.df, file=paste0(study,"_slides_all.df"))
     
#If you are going to load slides_all.df to save time, run in the command line:
#load(paste0(study,"_slides_all.df"))

### Make a spot annotations dataframe
annotation_targets.df <- filter(slides_all.df, slide_no==1, Block == 1 | Block == 2)
annotation_targets.df <- annotation_targets.df[,1:4]

annotation_targets.df <- cbind(row.names(annotation_targets.df), annotation_targets.df)
colnames(annotation_targets.df)[1] <- "target_id_numeric"
annotation_targets.df[6] <- c(paste(rownames(annotation_targets.df), annotation_targets.df$Name, annotation_targets.df$Block, sep = "_"))
colnames(annotation_targets.df)[6] <- "target_id_unique"
rownames(annotation_targets.df) <- c(annotation_targets.df$target_id_unique)

###Make a final index, this time of targets
index_target <- as.numeric(length(annotation_targets.df$target_id_numeric))

#####################################
###PROCESSING YOUR MICROARRAY data###
#####################################

###SELECT RELEVANT DATA###

###Generate specific data frames e.g. median foreground and background
#The column identifying the sample_id in slides_all.df id column 42
#The column identifying the foreground median in slides_all.df is 9
#The column identifying the background median in slides_all.df is 14
#The column identifying the median.df - the background in slides_all.df is 34 (here we ignore this column)

#Foreground 
if(iso == "IgM"){

fore.df <- annotation_targets.df
for(i in 1:length(samples)){
  ite_out<-slides_all.df$F488.Median[which(slides_all.df$Sample==samples_unique[i])]
  fore.df<-cbind(fore.df,ite_out)
  colnames(fore.df)[length(colnames(fore.df))]<-samples_unique[i]
}

#Background
back.df <- annotation_targets.df
for(i in 1:length(samples)){
  ite_out<-slides_all.df$B488.Median[which(slides_all.df$Sample==samples_unique[i])]
  back.df<-cbind(back.df,ite_out)
  colnames(back.df)[length(colnames(back.df))]<-samples_unique[i]
}

}

if(iso == "IgG"){
  
  fore.df <- annotation_targets.df
  for(i in 1:length(samples)){
    ite_out<-slides_all.df$F594.Median[which(slides_all.df$Sample==samples_unique[i])]
    fore.df<-cbind(fore.df,ite_out)
    colnames(fore.df)[length(colnames(fore.df))]<-samples_unique[i]
  }
  
  #Background
  back.df <- annotation_targets.df
  for(i in 1:length(samples)){
    ite_out<-slides_all.df$B594.Median[which(slides_all.df$Sample==samples_unique[i])]
    back.df<-cbind(back.df,ite_out)
    colnames(back.df)[length(colnames(back.df))]<-samples_unique[i]
  }
}


#Generate matrices of the same data (that is, the same data with only a single data type)
fore.matrix <- as.matrix(fore.df[,7:ncol(fore.df)])
back.matrix <- as.matrix(back.df[,7:ncol(back.df)])

###DATA CORRECTION###

### Perform background correction using limma function and normexp method.
cor.matrix <- backgroundCorrect.matrix(fore.matrix, back.matrix, method = "normexp", offset = 50, normexp.method = "mle")

#Export this for reference
write.csv(cor.matrix, file=paste0(study,"_background_corrected_MFI.csv"))

###Assign target names to groups of your array targets to identify their 'type'
targets_blank = c(grep("BLANK", annotation_targets.df$Name, ignore.case = TRUE))
targets_buffer = c(grep("buffer", annotation_targets.df$Name, ignore.case = TRUE))
targets_ref = c(grep("REF", annotation_targets.df$Name, ignore.case = TRUE))
targets_std = c(grep("Std", annotation_targets.df$Name, ignore.case = TRUE))
targets_allcontrol = c(targets_blank, targets_buffer, targets_ref, targets_std)

###GST subtraction!!! Do this with background corrected MFI before log transforming or anything.

#Prepare target data frame for merging (get "unique" names from annotation targets) 
target_meta2.df <- merge(target_meta.df, annotation_targets.df, by = "Name")

#merge target data frame with data
bunny <- merge(target_meta2.df, cor.matrix, by.y = "row.names", by.x = "target_id_unique", sort = FALSE, all.y = TRUE)
#subset based on GST tag

#depending on column name with tag upper or lowercase, use either of the next two lines
bun <- filter(bunny, Expression_tag == "GST")
#bun <- filter(bunny, Expression_Tag == "GST")
bun <- tibble::column_to_rownames(bun, var = "target_id_unique")

#isolate GST values and plot GST - 
GST <- bunny[grep("_GST_", bunny$target_id_unique),]

GSTmelt <- melt(GST)
GSTmelt <- filter(GSTmelt,!(variable == "Block" | variable == "Column" | variable == "Row"))

png(filename = paste0(study, "_GST_by_target_log.tif"), width = 7, height = 3.5, units = "in", res = 1200)

ggplot(GSTmelt, aes(x = variable, y=value, color = target_id_unique)) + geom_point() + theme_bw() +
  labs(x = "Sample", y = "Background Corrected MFI", title = "GST targets") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
  theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
  scale_y_continuous(breaks = 10**(1:10),trans = 'log10', labels = 10**(1:10))

graphics.off()

png(filename = paste0(study, "_GST.tif"), width = 5, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(2.1,4.1,2.1,2.1))
plot(c(GSTmelt$value), pch='*', col = "blue", main = "GST",
     ylab="Background Corrected MFI", xlab="GST spot", cex.main=1, cex.lab=1, cex.axis=0.7)

graphics.off()

if (reps == 1){
  GSTmean <- GST[,(ncol(target_meta2.df)+1):ncol(GST)]
}

if (reps == 2 | reps == 4){
  GSTmean <- colMeans(GST[,(ncol(target_meta2.df)+1):ncol(GST)])
}


#subtract GST directly for each sample from all tagged targets
#set negative values to 50, which is the minimum of the background corrected data (cor.matrix)

bundata <- bun[ncol(target_meta2.df):ncol(bun)]

subbedGST <- bundata
for(i in 1:ncol(subbedGST))
{
  subbedGST[,i] <- bundata[,i]-GSTmean[i]
  subbedGST[which(subbedGST[,i] <= 0),i] <- 50
}
remove(i)

#bind GST-subtracted data with data from non-tagged antigens, label final matrix as cor2.matrix.
#Make another data frame where the tagged protein values are replaced by their subtracted values
#filter out the GST tagged targets
no_tags.df <- filter(bunny, is.na(Expression_tag) | !(Expression_tag == "GST"))
#no_tags.df <- filter(bunny, is.na(Expression_Tag) | !(Expression_Tag == "GST"))
row.names(no_tags.df) <- no_tags.df$target_id_unique
no_tags.df <- no_tags.df[,(ncol(target_meta2.df)+1):ncol(no_tags.df)]

#then rbind the GST and the CD4 data frames to that one. 
cor1.matrix <- as.matrix(rbind(no_tags.df, subbedGST))

#need to put the matrix back in the same order as before (by target_id_unique) because of calling out buffer etc.
sortedcor <- merge(annotation_targets.df, cor1.matrix, by = "row.names", sort = FALSE)
cor2.matrix <- as.matrix(sortedcor[,8:ncol(sortedcor)])
row.names(cor2.matrix) <- row.names(annotation_targets.df)

#save ALL GST subtracted data in another file
write.csv(cor2.matrix, paste0(study, "_GST_subtracted_MFI.csv"))

###QUALITY CONTROL###

### 1. Background
#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples
#We are not automatically excluding samples with deviant background, but it might be a good idea to
# go back and review the images for those pads.

#Mean of median background MFI for all data points
back_mean <- mean(back.matrix)
#SD median background MFI fro all data points
back_sd <- sd(back.matrix)
#Cut-off for deviation from mean
back_cutoff <- back_mean+(3*back_sd)

###Identify deviant samples/slides/pads

#Generate a vector of mean background intensity for every target for EACH sample (background person magnitude)
back_sample_mean <- colMeans(back.matrix)
#Vectors to identify normal and deviant samples
back_normal <- which(back_sample_mean<=back_cutoff)
back_deviant <- which(back_sample_mean>back_cutoff)
#Generate a table showing which slides, pads, samples have deviant background values
back_sample_deviant <- samples.df[back_deviant,]
back_sample_deviant
write.csv(back_sample_deviant, file = paste0(study,"_deviant_sample_background.csv"))
#Coefficient of variation for all background datapoints, and all exlcuding deviant samples
back_cov_all <- sd(back.matrix)/mean(back.matrix)
back_cov_normal <- sd(back.matrix[, back_normal])/mean(back.matrix[, back_normal])

#Plots by slide/pad/sample
png(filename = paste0(study,"_background_mfi_samples.tif"), width = 5.5, height = 10, units = "in", res = 600)
par(mfrow=c(3,1), mar = c(2, 3, 2.25, 0.5), oma = c(11.5, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.6, cex.lab = 0.9, xpd=NA, las=1)
boxplot(t(back.matrix) ~ samples.df$slide_no, outcex=0.5,
        ylab="Background MFI", xlab="Slide", add=FALSE)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
title(main = "Median background MFI by SLIDE/SUBARRAY/SAMPLE\n", adj=0)
boxplot(t(back.matrix) ~ samples.df$block_rep_1, outcex=0.5,
        ylab="Background MFI", xlab="Subarrays (1/2, 3/4, etc.)", add=FALSE, las=1)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
boxplot(t(back.matrix) ~ samples.df$sample_id_unique, outcex=0.5,
        ylab="Background MFI", xlab="Samples", add=FALSE, las=2, cex.axis = 0.4, yaxt="n")
axis(2, cex.axis=0.6)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(c(paste("Mean overall background MFI:", round(back_mean, digits=3))), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("SD overall background MFI:", round(back_sd, digits=3))), side=1, cex=0.8, line=3.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Cut-off overall background MFI:", round(back_cutoff, digits=3))), side=1, cex=0.8, line=5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall background MFI:", round(back_cov_all, digits=3))), side=1, cex=0.8, line=6.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall background MFI (excl. deviant samples):", round(back_cov_normal, digits=3))), side=1, cex=0.8, line=8, outer=TRUE, xpd=NA, adj=0)

graphics.off()

###Identify deviant targets across all pads

#Mean background MFI for every sample for EACH target (background protein magnitude)
back_target_mean <- rowMeans(back.matrix)
#Mean background target magnitude for block 1, arranged in the order they are printed
back_target_mean_b1 = matrix(back_target_mean[annotation_targets.df$Block==1], nrow=max(annotation_targets.df$Row), ncol=max(annotation_targets.df$Column))
#Mean background target magnitude for block 2, arranged in the order they are printed
back_target_mean_b2 = matrix(back_target_mean[annotation_targets.df$Block==2], nrow=max(annotation_targets.df$Row), ncol=max(annotation_targets.df$Column))
#Are any background values for specific targets universally deviant, accross all pads?
back_target_deviant <- annotation_targets.df[back_target_mean>back_cutoff,]
back_target_deviant

#Plot background by target
png(filename = paste0(study,"_background_mfi_targets.tif"), width = 15, height = 10, units = "in", res = 600)
par(mfrow=c(2,1), mar = c(7, 3, 2.25, 0.5), oma = c(6, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.3, cex.lab = 0.9, xpd=NA, las=2)
boxplot(t(back.matrix), outcex=0.5, yaxt="n", ylab="Background MFI")
title(main="Background MFI by ARRAY TARGET: all samples", adj=0)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
axis(2, cex.axis=0.6)
mtext(c(paste("Targets 1-", index_target)), side=1, cex=0.9, xpd=NA, line=4.5, las=1)
boxplot(t(back.matrix[,back_normal]), outcex=0.5, yaxt="n", ylab="Background MFI")
title(main="Background median MFI by ARRAY TARGET: excl. deviant samples", adj=0)
axis(2, cex.axis=0.6)
mtext(c(paste("Targets 1-", index_target)), side=1, cex=0.9, xpd=NA, line=4.5, las=1)

mtext(c(paste("Mean overall background MFI:", round(back_mean, digits=3))), side=1, cex=0.8, line=0, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("SD overall background MFI:", round(back_sd, digits=3))), side=1, cex=0.8, line=1, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("Cut-off overall background MFI:", round(back_cutoff, digits=3))), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("CoV overall background MFI:", round(back_cov_all, digits=3))), side=1, cex=0.8, line=3, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("CoV overall background MFI (excl. deviant samples):", round(back_cov_normal, digits=3))), side=1, cex=0.8, line=4, outer=TRUE, xpd=NA, adj=0, las=1)

graphics.off()

#Cov sample (all samples and all targets)
png(filename = paste0(study,"_cov_back_sample.tif"), width = 10, height = 4, units = "in", res = 600)
par(mfrow=c(1,2), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

back_cov_sample <- c()
for(i in 1:ncol(back.matrix))
{
  temp <- sd(back.matrix[,i])/mean(back.matrix[,i])
  back_cov_sample <- c(back_cov_sample, temp)
}
plot(back_cov_sample, ylab="CoV background MFI", xlab="Sample", ylim=c(0,max(back_cov_sample)))
remove(temp)

back_cov_target <- c()
for(i in 1:nrow(back.matrix))
{
  temp <- sd(back.matrix[i,])/mean(back.matrix[i,])
  back_cov_target <- c(back_cov_target, temp)
}
plot(back_cov_target, ylab="CoV background MFI", xlab="Target", ylim=c(0,max(back_cov_target)))
remove(temp)
graphics.off()

### 2. Buffer
# **** Because we are no longer excluding buffer spots printed after certain targets, the different cor matrix used
#      and the different values for means, SD, etc are all the same. However, I have left both for ease of use 
#      because some values are used later in different plots and calculations

#Though some slides/pads may have high background, background correction may adjust this and make the data useable.
#The only way of identifying if this is true is by looking at the control spots by slide and pad.
#If the same samples have high controls, the background correction was not enough.
#Correction against plate buffers may be required...

#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The mean and cutoff will be calculated EXCLUDING any "bad" buffer targets previously set to NA.
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples

#***Samples which have deviant buffer means will be automatically excluded from further analysis
#***Buffer targets which are deviant across all pads and slides will be automatically excluded from further analysis

#Mean median corrected MFI for all data points
cor_buffer_mean <- mean(cor2.matrix[targets_buffer,], na.rm = TRUE)
#SD median corrected MFI for all data points
cor_buffer_sd <- sd(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Cut-off for deviation from mean
cor_cutoff <- cor_buffer_mean+(3*cor_buffer_sd)

###Identify deviant samples/slides/pads

# Generate a vector of mean buffer intensity for EACH sample (corrected person magnitude) 
cor_buffer_sample_mean <- colMeans(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Vectors to identify normal and deviant samples
cor_normal <- which(cor_buffer_sample_mean<=cor_cutoff)
cor_deviant <- which(cor_buffer_sample_mean>cor_cutoff)
#Generate a table showing which slides, pads, samples have deviant corrected buffer values
cor_sample_deviant <- samples.df[cor_deviant,]
cor_sample_deviant
write.csv(cor_sample_deviant, file = paste0(study,"_deviant_sample_buffer.csv"))

#Coefficient of variation for buffer datapoints and exlcuding deviant samples
cor_buffer_cov <- cor_buffer_sd/cor_buffer_mean
cor_buffer_cov_normal <- sd(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)/mean(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)

#Plot CoV of buffer by sample
cor_buffer_cov_sample <- c()

png(filename = paste0(study,"_cov_buffer.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

for(i in 1:ncol(cor2.matrix))
{
  temp <- sd(cor2.matrix[targets_buffer,i], na.rm = TRUE)/mean(cor2.matrix[targets_buffer,i], na.rm = TRUE)
  cor_buffer_cov_sample <- c(cor_buffer_cov_sample, temp)
}
plot(cor_buffer_cov_sample, ylab="CoV corrected buffer MFI", xlab="Sample", ylim=c(0,max(cor_buffer_cov_sample)))
remove(temp)
graphics.off()

#Buffer assessment INCLUDING "bad" spots:
cor_buffer_all_mean <- mean(cor.matrix[targets_buffer,])
#SD median corrected MFI for all data points
cor_buffer_all_sd <- sd(cor.matrix[targets_buffer,])
#Cut-off for deviation from mean
cor_all_cutoff <- cor_buffer_all_mean+(3*cor_buffer_all_sd)

#Coefficient of variation for all buffer datapoints (INCLUDING "bad" spots), and all exlcuding deviant samples
cor_buffer_cov_all <- cor_buffer_all_sd/cor_buffer_all_mean
cor_buffer_cov_all_normal <- sd(cor.matrix[targets_buffer, cor_normal])/mean(cor.matrix[targets_buffer, cor_normal])

#Plots by slide/pad/sample INCLUDING "bad" spots (all buffer data)
png(filename = paste0(study,"_buffer_mfi_QCplots.tif"), width = 5.5, height = 10, units = "in", res = 600)
par(mfrow=c(3,1), mar = c(2, 4, 2.25, 0.5), oma = c(11.5, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 1, cex.lab = 1.25, xpd=NA, las=1)

boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$slide_no, outcex=0.5, xlab="Slide", add=FALSE, log = "y")
abline(h = cor_all_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
title(main = "Corrected Buffer MFI by SLIDE/SUBARRAY/SAMPLE\n", adj=0)
title(ylab="Corrected MFI (log scale)", line=2.7)

boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$block_rep_1, outcex=0.5, xlab="Subarrays (1/2, 3/4, etc.)", add=FALSE, las=1, log = "y")
abline(h = cor_all_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
title(ylab="Corrected MFI (log scale)", line=2.7)

boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$sample_id_unique, outcex=0.5, xlab="Sample", add=FALSE, las=2, cex.axis = 0.4, yaxt="n", log = "y")
axis(2, cex.axis=1)
abline(h = cor_all_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
abline(h = cor_buffer_all_mean, col = "red", lwd = 0.7, xpd=FALSE)
title(ylab="Corrected MFI (log scale)", line=2.7)

mtext(c(paste("Mean overall corrected buffer MFI:", round(cor_buffer_all_mean, digits=3))), side=1, cex=0.8, line=3.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("SD overall corrected buffer MFI:", round(cor_buffer_all_sd, digits=3))), side=1, cex=0.8, line=5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Cut-off overall corrected buffer MFI:", round(cor_all_cutoff, digits=3))), side=1, cex=0.8, line=6.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI:", round(cor_buffer_cov_all, digits=3))), side=1, cex=0.8, line=8, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI (excl. deviant samples):", round(cor_buffer_cov_all_normal, digits=3))), side=1, cex=0.8, line=9.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

###Identify deviant buffer targets across all pads

#Mean corrected MFI for every sample for EACH target (background protein magnitude)
cor_target_mean <- rowMeans(cor2.matrix)
#Are any corrected buffer targets universally deviant, accross all pads?
deviant_buffer_targets <- Reduce(intersect, list(targets_buffer, which(cor_target_mean>cor_cutoff)))
cor_buffer_deviant.df <- annotation_targets.df[deviant_buffer_targets,]
cor_buffer_deviant.df
write.csv(cor_buffer_deviant.df, file = paste0(study,"_deviant_buffer_targets.csv"))

#list of buffer targets removed because they were deviant (mean > cutoff)
removed_buffer_targets <- c()
for (i in 1:length(targets_buffer)){
  for(j in 1:nrow(cor_buffer_deviant.df))
    if (targets_buffer[i] == cor_buffer_deviant.df[j,1]){
      removed_buffer_targets[j] <- targets_buffer[i]
    }
}
remove(i,j)
removed_buffer_targets <- subset(removed_buffer_targets, !is.na(removed_buffer_targets))

#All slides, INCLUDING "bad" spots (ALL buffer targets)
png(filename = paste0(study, "_buffer_targets.tif"), width = 11, height = 4, units = "in", res = 600)
par(mar = c(5, 3, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 1, xpd=NA, las=2)
boxplot(t(cor.matrix[targets_buffer,]),
        cex=0.5,
        ylab="Corrected MFI (log scale)", log = "y")
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(paste("Deviant buffer targets:", paste(removed_buffer_targets, collapse = ",")), las = 1)

graphics.off()

#By slide for up to 6 slides (the last slides) - INCLUDING "bad" spots (ALL buffer targets)
#mfrow is ordered by the number of rows, and columns of plots you want - so must be edited based on the number of slides
png(filename = paste0(study, "_buffer_spots_slide.tif"), width = 5, height = 4, units = "in", res = 600)
par(mfrow=c(2,3), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=2)
for (i in c(1, 20, 40, 60, 80, 98)){
  boxplot(t(cor.matrix[targets_buffer,samples.df$slide_no==i]),
          ylab="Corrected MFI",
          ylim=c(0,2000),
          add=FALSE, 
          cex=0.5,
          xpd=NA,
          main=c(paste("Slide", i)),
          adj=0)
  abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
}

graphics.off()

#Automatically set to NA the samples with deviant buffer values in a new matrix
cor3.matrix <- cor2.matrix
cor3.matrix[,cor_deviant] <- NA

#Add setting exclude = yes for these samples in the sample meta data frame
for(i in 1:nrow(cor_sample_deviant)){
  sample_meta.df$exclude[which(sample_meta.df$sample_id == cor_sample_deviant$sample_id[i])] <- "yes"
}

#Automatically set to NA the buffer targets deviant across all arrays
# Do this even though we are also removing outliers because this 
# indicates a systematic issue with that spot which means it's data should 
# not be counted because it could skew the mean used for buffer normalization
cor3.matrix[deviant_buffer_targets,] <- NA

#for Keneba study, remove the 4 systematically different buffer spots
#due to unknown problem with printing after menX
sysdeviants <- c("102_buffer_1", "114_buffer_1", "390_buffer_2", "402_buffer_2")
cor3.matrix[sysdeviants,] <- NA

### Mean values for each target ordered by position within the arrays (not log-transformed or normalized data)
#KG - I think we might want this somehwere else in the script using a more processed matrix

#Check average corrected values for each target for all individuals
cor_target_mean <- rowMeans(cor3.matrix, na.rm = TRUE)
#Mean background target magnitude for block 1, arranged in the order they are printed
cor_target_mean_b1 = t(matrix(round(cor_target_mean[annotation_targets.df$Block==1], digits=2), nrow=max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
#Mean background target magnitude for block 2, arranged in the order they are printed
cor_target_mean_b2 = t(matrix(round(cor_target_mean[annotation_targets.df$Block==2], digits=2), nrow=max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
cor_target_mean_b1b2 <- rbind(cor_target_mean_b1, cor_target_mean_b2)
#Annotation plate maps
annotation_target_b1 <- t(matrix(annotation_targets.df$target_id_unique [annotation_targets.df$Block==1], nrow = max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
annotation_target_b2 <- t(matrix(annotation_targets.df$target_id_unique [annotation_targets.df$Block==2], nrow = max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
annotation_target_b1b2 <- rbind(annotation_target_b1, annotation_target_b2)
#Write csv, which can be presented as a heatmap
write.csv(cbind(cor_target_mean_b1b2, annotation_target_b1b2), file=paste0(study,"_target_mean_as_array.csv"))
remove(cor_target_mean, cor_target_mean_b1, cor_target_mean_b2,cor_target_mean_b1b2, annotation_target_b1, annotation_target_b2, annotation_target_b1b2)

### Remove Buffer Outliers - this section copied from Sanger v2 data processing.R in the master branch.

#plot all buffer values as a histogram or qplot to check normality - only normal after log2 transformation, not MFI
png(filename = paste0(study, "_Buffer_hist.tif"), width = 5, height = 7.5, units = "in", res = 1200)
par(mfrow=c(2,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

hist(c(cor3.matrix[targets_buffer,]), breaks = 25, col = "blue",
     ylab="Frequency", xlab="Background corrected MFI", main = NULL)
title(main="All Buffer, MFI", adj=0)

hist(log2(c(cor3.matrix[targets_buffer,])), col = "darkblue", breaks = 25,
     ylab="Frequency", xlab="Log2(MFI)", main = NULL)
title(main="All Buffer, Log2(MFI)", adj=0)

graphics.off()

#isolate and transform buffer values
Buffer <- cor3.matrix[targets_buffer,]
Buflog <- log2(Buffer)

#check the 1.5*IQR method instead? 
#set rule = 1.5 initially - this was getting 6.15, way too low
#The problem is the IQR is only 0.22 because the data is mostly all at the bottom end
#So this method will not work. 
maxIQR <- function(data, rule){
  x <- quantile(data,0.75, na.rm=TRUE) + (IQR(data, na.rm=TRUE) * rule)
  return(x)
}

outliers <- maxIQR(Buflog, 3) #this is getting a result of a cutoff of 8.385736 on the log2 data (Keneba pilot IgG) and 6.924607 IgM, which is too low

#Altered the scores function source code so that it will be able to deal with NAs 
#verified to work for the z type, hypothetically for the other types as well.

"scores" <- function (x, type = c("z","t","chisq","iqr","mad"), prob = NA, lim = NA) 
  {
    if (is.matrix(x)) 
      apply(x, 2, scores, type = type, prob = prob, lim = lim)
    else if (is.data.frame(x)) 
      as.data.frame(sapply(x, scores, type = type, prob = prob, lim = lim))
    else {
      n <- length(x)
      s <- match.arg(type)
      ty <- switch(s, z=0, t=1, chisq=2, iqr=3, mad=4)
      
      if (ty == 0) {
        res <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
        if (is.na(prob)) res
        else {
          if (prob == 1) pnorm(res)
          else	if (prob == 0) abs(res) > (n-1)/sqrt(n)
          else abs(res) > qnorm(prob)
        }
      }
      else if (ty == 1) {
        t <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
        res <- (t*sqrt(n-2))/sqrt(n-1-t^2)
        if (is.na(prob)) res
        else {
          if (prob == 1) pt(res,n-2)
          else	if (prob == 0) abs(res) > (n-1)/sqrt(n)
          else abs(res) > qt(prob,n-2)
        }
        
      }
      else if (ty == 2) {
        res <- (x - mean(x, na.rm = TRUE))^2/var(x, na.rm = TRUE)
        if (is.na(prob)) res
        else {
          if (prob == 1) pchisq(res,1)
          else abs(res) > qchisq(prob,1)
        }
      }
      else if (ty == 3) {
        res <- x
        Q1 <- quantile(x,0.25, na.rm = TRUE)
        Q3 <- quantile(x,0.75, na.rm = TRUE)
        res[x >= Q1 & res <= Q3] <- 0
        res[x < Q1] <- (res[x < Q1]-Q1)/IQR(x)
        res[x > Q3] <- (res[x > Q3]-Q3)/IQR(x)
        if (is.na(lim)) res
        else abs(res) > lim
      }
      else if (ty == 4) {
        res <- (x - median(x, na.rm = TRUE))/mad(x, na.rm = TRUE)
        if (is.na(prob)) res
        else {
          if (prob == 1) pnorm(res)
          else	if (prob == 0) abs(res) > (n-1)/sqrt(n)
          else abs(res) > qnorm(prob)
        }
      }
      
    }
  }

#calculate outliers for log2 transformed data, for all buffer data considered as one population
#decided to go with prob = 0.99 for Keneba Pilot 
#although this is a conservative cutoff based on the QC plots
BUFout <- scores(Buflog, type = "z", prob = 0.99)

#remove outliers (set to NA in original background corrected MFI data) 
BUFrm995 <- as.matrix(Buffer)
for(i in 1:length(BUFrm995)){
  if(is.na(BUFout[[i]])){ 
    i = i+1 
  }else if(BUFout[[i]] == TRUE){
    BUFrm995[[i]] <- NA
  }
}

max(BUFrm995, na.rm = TRUE) #843 for 0.99 (IgG main study v1), 401.1728 (IgG main study v2)
#117 for 0.99 for IgM
min(BUFrm995, na.rm = TRUE) #50 for 0.99 (IgG main study v1 and v2)
#50 for 0.99 for IgM

#what percentage of buffer values are removed?
length(which(BUFout == TRUE))/length(BUFout)*100 # 2.439413% for 0.99 IgG main study v1
#2.063138 for IgG main study v2
#2.976941% for 0.99 for IgM 

#plot histogram again with outliers removed 
png(filename = paste0(study, "_Buffer_hist_p.99.tif"), width = 5, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

hist(log2(c(BUFrm995)), breaks = 25, col = "blue",
     ylab="Frequency", xlab="Log2(MFI)", main = NULL)
title(main="Buffer without outliers, MFI", adj=0)

graphics.off()

#in the main cor matrix, put in buffer matrix with outliers removed
nobuffer <- cor3.matrix[-targets_buffer,]
newbuffers <- rbind(nobuffer, BUFrm995)

#need to put the matrix back in the same order as before (by target_id_unique) because of calling out buffer etc.
sortedcorbuf <- merge(annotation_targets.df, newbuffers, by = "row.names", sort = FALSE)
cor4.matrix <- as.matrix(sortedcorbuf[,8:ncol(sortedcorbuf)])
row.names(cor4.matrix) <- row.names(annotation_targets.df)

#### LOG TRANSFORMATION AND NORMALIZATION ###

### Log transform the data (base 2)
log.cor.matrix <- log2(cor4.matrix)

#export this matrix to compare normalized vs. not normalized data.
write.csv(t(log.cor.matrix), file = paste0(study,"_log_data.csv"))

### Normalization

###Create sample specific buffer means for normalisation
cor2_buffer_sample_mean <- colMeans(cor4.matrix[targets_buffer,], na.rm = TRUE)
cor2_buffer_sample_sd <- apply(cor4.matrix[targets_buffer,], 2, sd, na.rm = TRUE)
#log2 transform sample buffer mean
log_buffer_sample_mean <- log2(cor2_buffer_sample_mean)

#subtract buffer mean from each sample to generate new intensity matrices (log transformed data)
norm.matrix <- log.cor.matrix
for(i in 1:ncol(norm.matrix))
{
  norm.matrix[,i] <- norm.matrix[,i]-log_buffer_sample_mean[i]
}

write.csv(t(norm.matrix), file = paste0(study,"_normalized_log_data.csv"))

### PLOTTING STANDARDS ###

### Plot standard values for each sample and assess variation - with negative normalized values

  #isolate data for standards, normalized and not normalized
  stds_norm <- norm.matrix[targets_std,]
  stds_pre <- log.cor.matrix[targets_std,]
  
  #isolate data from isotype probed
  isostd_norm <- stds_norm[grep(iso, row.names(stds_norm)),]
  isostd_pre <- stds_pre[grep(iso, row.names(stds_pre)),]

  #Plot Std 3 in Levey Jennings Style plot
  #KG - I don't know if these column numbers will hold up every time...switch to using grep?
  #average replicates for reps == 2, arithmetic mean!
  if (reps == 1){
    std_3_norm <- isostd_norm[c(3),]
    std_3_pre <- isostd_pre[c(3),]
    }

  if (reps == 2){
    std_3_norm <- isostd_norm[c(3,9),]
    std_3_norm <- log2(apply((2^std_3_norm), 2, mean))
  
    std_3_pre <- isostd_pre[c(3,9),]
    std_3_pre <- log2(apply((2^std_3_pre), 2, mean))
  }
  
  if (reps == 4){
    std_3_norm <- isostd_norm[c(3,9,15,21),]
    std_3_norm <- log2(apply((2^std_3_norm), 2, mean))
    
    std_3_pre <- isostd_pre[c(3,9,15,21),]
    std_3_pre <- log2(apply((2^std_3_pre), 2, mean))
  }

  #calculate geometric (not arithmetic) mean, SD and CV
  #normalized:
  std3mean <- mean(c(std_3_norm), na.rm = TRUE)
  std3sd <- sd(c(std_3_norm), na.rm = TRUE)
  e_std3sd <- std3sd*log(2)
  std3cv <- sqrt(exp(e_std3sd^2)-1)*100

  #pre-normalized:
  std3mean1 <- mean(c(std_3_pre), na.rm = TRUE)
  std3sd1 <-sd(c(std_3_pre), na.rm = TRUE)
  e_std3sd1 <- std3sd1*log(2)
  std3cv1 <- sqrt(exp(e_std3sd1^2)-1)*100

  #Plotting Std 3 Levey Jennings Style
  png(filename = paste0(study, "_std_3_LJ.tif"), width = 5, height = 7.5, units = "in", res = 1200)
  par(mfrow=c(2,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  plot(c(std_3_norm), pch='*', col = "blue", ylim=c(min(std_3_norm, na.rm = TRUE),max(std_3_norm, na.rm=TRUE)),
     ylab="Normalized log2(MFI)", xlab="Sample (Array)")

  abline(h=std3mean)
  abline(h=std3mean+2*std3sd,lty=2)
  abline(h=std3mean-2*std3sd,lty=2)
  abline(h=std3mean+std3sd,lty=3)
  abline(h=std3mean-std3sd,lty=3) 

  plot(c(std_3_pre), pch='*', col = "darkblue", ylim=c(min(std_3_pre, na.rm = TRUE),max(std_3_pre, na.rm=TRUE)),
     ylab="log2(MFI) (NOT normalized)", xlab="Sample (Array)")

  abline(h=std3mean1)
  abline(h=std3mean1+2*std3sd1,lty=2)
  abline(h=std3mean1-2*std3sd1,lty=2)
  abline(h=std3mean1+std3sd1,lty=3)
  abline(h=std3mean1-std3sd1,lty=3)

  mtext(paste("Geometric CV, Normalized:", round(std3cv, digits=2), "%" ), side=1, cex=0.8, line=0.5, outer=TRUE, xpd=NA, adj=0)
  mtext(paste("Geometric CV, NOT Normalized:", round(std3cv1, digits=2), "%"), side=1, cex=0.8, line=1.5, outer=TRUE, xpd=NA, adj=0)

  graphics.off()

###Plot All standards together in ggplot2, but separately for rep1 and rep2 - for ALL THREE std curves
  
  Ig <- c("IgG Std", "IgM Std", "IgM Fc")
  
for(i in 1:length(Ig)){
    
  type = Ig[i]
  
  norm <- stds_norm[grep(type, row.names(stds_norm)),]
  pre <- stds_pre[grep(type, row.names(stds_pre)),]
  
  #normalized
  std_norm_1 <- norm[1:6,]
  std_norm_2 <- norm[7:12,]
  
  std1melt <- melt(std_norm_1, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_norm_1_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std1melt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep1 Normalized") + ylim(-1,10) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
  std2melt <- melt(std_norm_2, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_norm_2_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std2melt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep2 Normalized") + ylim (-1,10) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
  #not normalized
  std_pre_1 <- pre[1:6,]
  std_pre_2 <- pre[7:12,]
  
  std1premelt <- melt(std_pre_1, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_pre_1_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std1premelt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep1 Pre-Normalization") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
  std2premelt <- melt(std_pre_2, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_pre_2_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std2premelt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep2 Pre-Normalization") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
  }
  
#Do this plot again with one plot each isotype, for the mean of all the standard reps
#TBD

  
### Set negative normalized log values to zero. This will be used for some analyses. 
# For other analyses, including sending data to Nuno, the data will be input without setting values to 0. 
# From now on, the norm.matrix has negative values, and norm2.matrix does not. 
norm2.matrix <- norm.matrix

for (i in 1:length(norm2.matrix))
{
  if (is.na(norm2.matrix[[i]]) | norm2.matrix[[i]] > 0) 
  {
    i = i+1
  } else if(norm2.matrix[[i]] < 0) 
  {
    norm2.matrix[[i]] <- 0
    i = i+1
  }
}
remove(i)

write.csv(t(norm2.matrix), file = paste0(study,"_normalized_log_data_0s.csv"))


### Average duplicates - INCLUDING negative values, if the data has technical replicates in the form of 2 blocks / subarray

if (reps == 2)
{
  n = nrow(norm.matrix)/2
  rep1 <- norm.matrix[1:n,]
  rep2 <- norm.matrix[(n+1):(n*2),]
  
  normaverage.matrix <- matrix(nrow = n, ncol = ncol(norm.matrix))
  colnames(normaverage.matrix) = colnames(norm.matrix)
  rownames(normaverage.matrix) = rownames(norm.matrix[(1:n),])
  
  normaverage.matrix <- log2((2^rep1 + 2^rep2)/2)
  
}

### Average quadruplicates
if (reps == 4)
{
  n = nrow(norm.matrix)/2
  block1 <- norm.matrix[1:n,]
  block2 <- norm.matrix[(n+1):(n*2),]
  
  #calling both the reps in block 1 reps 1 and 2 (144 total spots per replicate)
  #the first 12 are rep 1, second 12 rep 2, etc so make a number sequence
  rep1num <- c(1, (1+ cumsum(rep(c(1,1,1,1,1,1,1,1,1,1,1,13), 12))))
  rep1 <- block1[rep1num[1:144],]
  
  rep2num <- c(13, (13 + cumsum(rep(c(1,1,1,1,1,1,1,1,1,1,1,13), 12))))
  rep2 <- block1[rep2num[1:144],]
  
  #the reps in block 2 are reps 3 and 4
  rep3 <- block2[rep1num[1:144],]
  rep4 <- block2[rep2num[1:144],]
  
  normaverage.matrix <- matrix(nrow = nrow(rep1), ncol = ncol(norm.matrix))
  colnames(normaverage.matrix) = colnames(norm.matrix)
  rownames(normaverage.matrix) = rownames(rep1)
  
  normaverage.matrix <- log2((2^rep1 + 2^rep2 + 2^rep3 + 2^rep4)/4)
  
}


### Check for deviant technical replicates, automatically exclude (set to NA)
# Use a modified Patrick's formula for ELISA called Katie's formula for microarray ;) to compare replicates within one array
  # if rep1 or rep2 is more than 2 times rep2 or rep1, respectively, exclude that pair
  # only apply the test if the values for both rep1 and rep2 are above 2 (on log2 scale)
  # Also, can redo this using the subsetted matrices (rep1 and rep2) and it should be shorter
if (reps == 2)
{ 
  for (k in 1:ncol(rep1))
  {
    for(j in 1:nrow(rep1)) 
    {
      if(is.na(rep1[j,k]) | is.na(rep2[j,k]) | (rep1[j,k]<2 & rep2[j,k]<2)){
        j+1
      } else if (rep1[j,k] > (log2(2) + rep2[j,k]) | (rep2[j,k] > (log2(2) + rep1[j,k])) == TRUE) 
      {
        normaverage.matrix[j,k] <- NA
      }
    }
  }
  remove(j,k)
  
  write.csv(normaverage.matrix, paste0(study, "_average_norm_log_data.csv")) 
  
  ## Calculate correlation coefficient (default is pearson). Deviants are still included.
  repR <- cor(c(rep1), c(rep2), use = "complete.obs")
  repR34 <- cor(c(rep3), c(rep4), use = "complete.obs")
  
  print(repR)
  
  ## Plot replicate 1 v. replicate 2 for each protein or each person and calculate correlation coefficient.
  png(filename = paste0(study, "_replicatescorrelation.tif"), width = 5, height = 4, units = "in", res = 600)
  par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
      mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
  
  plot(rep1, rep2, col="red", cex = 0.1)
  mtext(c(paste("Pearson correlation coefficient:", round(repR, digits=4))), side=3, adj=0)
  
  graphics.off()
  
  ## Plot replicate 3 v. replicate 4 for each protein or each person and calculate correlation coefficient.
  png(filename = paste0(study, "_replicatescorrelation34.tif"), width = 5, height = 4, units = "in", res = 600)
  par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
      mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
  
  plot(rep3, rep4, col="red", cex = 0.1)
  mtext(c(paste("Pearson correlation coefficient:", round(repR34, digits=4))), side=3, adj=0)
  
  graphics.off()
  
}

### correlation plot of all the reps vs each other
#make a matrix of all the data

if (reps == 4){
  repsmat <- cbind(c(rep1), c(rep2), c(rep3), c(rep4))
  colnames(repsmat) <- c("Rep 1", "Rep 2", "Rep 3", "Rep 4")
  
  png(filename = paste0(study, "_Rep_Correlogram.tif"), width = 7, height = 6.5, units = "in", res = 1200)
  
  corrplot.mixed(cor(repsmat, use = "complete.obs"), tl.col="black", tl.pos = "lt")
  
  graphics.off()
  
}


### Average duplicates - Negative values set to 0s, if the data has technical replicates in the form of 2 blocks / subarray

if (reps == 2)
{
  n = nrow(norm2.matrix)/2
  rep1 <- norm2.matrix[1:n,]
  rep2 <- norm2.matrix[(n+1):(n*2),]
  
  norm2average.matrix <- matrix(nrow = n, ncol = ncol(norm2.matrix))
  colnames(norm2average.matrix) = colnames(norm2.matrix)
  rownames(norm2average.matrix) = rownames(norm2.matrix[(1:n),])
  
  norm2average.matrix <- log2((2^rep1 + 2^rep2)/2)
  
}

### Check for deviant technical replicates, automatically exclude (set to NA)
# Use Patrick’s formula for ELISA to compare replicates within one array
# if rep1 or rep2 is more than 1.5 times rep2 or rep1, respectively, exclude that pair
# Also, can redo this using the subsetted matrices (rep1 and rep2) and it should be shorter
if (reps == 2)
{ 
  for (k in 1:ncol(rep1))
  {
    for(j in 1:nrow(rep1)) 
    {
      if(is.na(rep1[j,k]) | is.na(rep2[j,k]) | (rep1[j,k]<2 & rep2[j,k]<2)){
        j+1
      } else if (rep1[j,k] > (log2(2) + rep2[j,k]) | (rep2[j,k] > (log2(2) + rep1[j,k])) == TRUE) 
      {
        norm2average.matrix[j,k] <- NA
      }
    }
  }
  remove(j,k)
  
  write.csv(norm2average.matrix, paste0(study, "_average_norm_log_data_0s.csv")) 
  
  ## Calculate correlation coefficient (default is pearson). Deviants are still included.
  repR <- cor(c(rep1), c(rep2), use = "complete.obs")
  print(repR)
  
  ## Plot replicate 1 v. replicate 2 for each protein or each person and calculate correlation coefficient.
  png(filename = paste0(study, "_replicatescorrelation_0s.tif"), width = 5, height = 4, units = "in", res = 600)
  par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
      mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
  
  plot(rep1, rep2, col="red", cex = 0.1)
  mtext(c(paste("Pearson correlation coefficient:", round(repR, digits=4))), side=3, adj=0)
  
  graphics.off()
  
}

#create matrix that is the same name whether or not we needed to average duplicates
#Including negative values
if (reps==2){norm3.matrix<-normaverage.matrix}
if (reps==1){norm3.matrix<-norm.matrix}

#With negative values set to 0
if (reps==2){norm4.matrix<-norm2average.matrix}
if (reps==1){norm4.matrix<-norm2.matrix}

###Identifying and excluding samples assayed in duplicate on the arrays 

#Create a transposed version of the data matrix (includes negative values)
trans.norm.matrix <- t(norm3.matrix)
for_dups.matrix <- tibble::rownames_to_column(as.data.frame(trans.norm.matrix), var = "sample_id_unique")

#Identify duplicate samples and store results in a data frame
#This includes the "test" sample type only, not the controls!
dup_samples <- samples.df[(duplicated(samples.df$sample_id) | duplicated(samples.df$sample_id, fromLast=TRUE)),]
dup_samples <- dup_samples[dup_samples$sample_type=="test",]

dup_samp_data <- merge(dup_samples, for_dups.matrix, by = "sample_id_unique")

#Export a table with all the info and data for the duplicates only
write.csv(dup_samp_data, file = paste0(study, "_duplicate_assayed_samples.csv"))

#Set exclude to "yes" for all samples assayed in duplicate.
#for(i in 1:nrow(dup_samples)){
#  sample_meta.df$exclude[which(sample_meta.df$sample_id == dup_samples$sample_id[i])] <- "yes"
#}

###Exporting processed data and metadata for further analysis (i.e. to give to Nuno)
#This data includes negative normalized values.

#Character vector of samples to be removed
samples_exclude <- sample_meta.df$sample_id_unique[which(sample_meta.df$exclude =="yes")]

#Sample metadata file, with samples removed if exclude == yes. 
  #Sample metadata and normalized log data are linked by the column "sample_id_unique"
  sample_meta_f.df <- sample_meta.df[(!(sample_meta.df$exclude == "yes") | is.na(sample_meta.df$exclude)),]

  #Export file
  write.csv(sample_meta_f.df, file = paste0(study, "_sample_metadata.csv"))

#Normalized log data with samples as rows and targets as columns for every target(including controls)

  #With samples_exclude removed and convert to data frame
  trans.norm2.matrix <- trans.norm.matrix[(!rownames(trans.norm.matrix) %in% samples_exclude),]
  trans.norm.df <- as.data.frame(trans.norm2.matrix)
  
  #Change the rownames to a separate column - 1st column is "sample_id_unique"
  trans.norm.df <- tibble::rownames_to_column(trans.norm.df, var = "sample_id_unique")
  
  #Export file
  write.csv(trans.norm.df, file = paste0(study, "_final_processed_data.csv"))

#Target metadata for every target - nothing has changed about this since the beginning
  #Export file
  write.csv(target_meta.df, file = paste0(study, "_target_metadata.csv"))

#Prepare final data with GST subtracted, control targets removed - For Negs set to 0 ONLY! (norm4.matrix)
  #Assign sample type names, to identify control and test samples (logical)
  samples_test <- sample_meta.df$sample_id_unique[which(sample_meta.df$sample_type =="test")]
  samples_control <- sample_meta.df$sample_id_unique[which(sample_meta.df$sample_type =="control")]
  
  #Define a list of targets to be removed from further analysis (controls)
  rmsamp_all <- unique(c(targets_blank, targets_buffer, targets_ref, targets_std))
  
  #Remove control protein targets - Don't remove control samples yet, need to do tag subtraction 
  #from those samples as well, and want them included in some exported data.
  #Do remove samples that should be excluded
  norm_sub.matrix <- norm4.matrix[-rmsamp_all,(!colnames(norm4.matrix) %in% samples_exclude)]
  
  #Replace current target names with original target names now that control targets are removed
  #might be useful to merge this instead with the target dataframe?
  norm_sub3.df <- merge(norm_sub.matrix, annotation_targets.df, by ="row.names", sort = FALSE)
  norm_sub3.df <- tibble::column_to_rownames(norm_sub3.df, var="Row.names")
  row.names(norm_sub3.df) <- norm_sub3.df$Name
  norm_sub4.df <- norm_sub3.df[,1:ncol(norm_sub.matrix)]
  
  #Make the dilution column of target_meta.df a character type
  #target_meta.df$Concentration <- as.character(target_meta.df$Concentration)
  
  #Merge with target metadata to filter based on expression tag etc.
  target.df <- merge(target_meta2.df, norm_sub4.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
  
  #Save final data frame to csv file
  write.csv(norm_sub4.df, paste0(study, "_finalafterprocessing.csv"))
  
#Save R workspace so that can load prior to analysis 
  save.image(file= paste0(study,"_AfterProcessing.RData"))

  
  
  