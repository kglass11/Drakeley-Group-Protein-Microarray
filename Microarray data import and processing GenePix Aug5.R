###Combined script for reading in and processing microarray data to prepare for analysis
  #IgG or IgM
  #Last updated August 5, 2018

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

### Define variables based on your study that will be used later in the script
# define working directory character vector, example "I:/Drakeley Group/Protein microarrays/Experiments/030417 Ghanaian samples/RepeatProcessingMay21KG"
workdir <- "I:/Drakeley Group/PROTEIN MICROARRAYS/Experiments/Elin MSC Project"

# define a shorthand name for your study which will be appended in the file name of all exported files
study <- "CHMI IgG"

#define file name for sample IDs character vector, example "Analysis sample list 2.csv"
sample_file <- "CHMI_microarray_samples.csv"

#define file name for sample file + additional metadata (character vector)
meta_file <- "Sample_metadata.csv"

#define file name for antigen list file with additional info about targets.
target_file <- "Target metadata.csv" 

#number of technical replicates for the study (usually 1 or 2)
reps <- 2

#define number of blocks per slide
index_block <- 32

#define isotype
iso <- "IgG"

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

###Bind all data from the slide data list (slides.list) into a single dataframe
#you may get a warning after this step, invalid factor level, this is not a problem!
slides_all.df <- c()
for(i in 1:length(slides_list)) { 
  
  slides_all.df <- rbind(slides_all.df,slides_list[[i]])
}
remove(i)

###Read in list of sample IDs
samples.df <- read.csv(sample_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

###Read in sample metadata file
sample_meta1.df <- read.csv(meta_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

###Read in target metadata file
target_meta.df <- read.csv(target_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

###Processing sample list:
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

#Export a table of duplicate entries that do not match (i.e. there is a problem with epi data) 
duplicate_metadata <- sample_meta2.df[(duplicated(sample_meta2.df$sample_id)| duplicated(sample_meta2.df$sample_id, fromLast=TRUE)),]
write.csv(duplicate_metadata, file = paste0(study, "_duplicate_metadata.csv"))

#Remove duplicate entries automatically (both duplicates)
sample_meta3.df <- sample_meta2.df[!(duplicated(sample_meta2.df$sample_id) | duplicated(sample_meta2.df$sample_id, fromLast=TRUE)),]

#Set exclude = "yes" in sample list file for duplicates that don't match in metadata
dup <- unique(duplicate_metadata$sample_id)
for(i in 1:length(dup)){
  samples.df$exclude[which(samples.df$sample_id == dup[i])] <- "yes"
}
remove(i)

### Merge the sample list file and the sample metadata file to include the appropriate metadata
#The duplicate metadata will now be listed as NA, with exclude = yes

sample_meta.df <- merge(samples.df, sample_meta3.df, by = "sample_id", all.x = TRUE)

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

#*** Old code for the same thing:
###Assign your sample_ids to each row of the combined slide data (slides_all.df)
#The order of data in your samples.df file is irrelevant, as long as each sample ID is correctly matched to its slide and block numbers

# for(i in 1:dim(slides_all.df)[1]){
#   
#   print(i)
#   row_ite<-slides_all.df[i,]
#   block_ite<-row_ite$Block
#   slide_ite<-row_ite$slide_no
#   sample_info_1<-samples.df[which(samples.df$slide_no==slide_ite),]
#   match<-block_ite%in%sample_info_1[,"block_rep_1"]
#   if(match==TRUE){
#     value<-which(block_ite==sample_info_1[,"block_rep_1"])
#   }else{
#     value<-which(block_ite==sample_info_1[,"block_rep_2"])
#   }
#   sample_info_2<-sample_info_1[value,"sample_id_unique"]
#   slides_all.df$Sample[i]<-as.character(sample_info_2)
# }
# remove(i, sample_info_1, sample_info_2, match, row_ite, block_ite, slide_ite)

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
if(iso == "IgG"){

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

if(iso == "IgM"){
  
  fore.df <- annotation_targets.df
  for(i in 1:length(samples)){
    ite_out<-slides_all.df$F635.Median[which(slides_all.df$Sample==samples_unique[i])]
    fore.df<-cbind(fore.df,ite_out)
    colnames(fore.df)[length(colnames(fore.df))]<-samples_unique[i]
  }
  
  #Background
  back.df <- annotation_targets.df
  for(i in 1:length(samples)){
    ite_out<-slides_all.df$B635.Median[which(slides_all.df$Sample==samples_unique[i])]
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

###Remove "bad" spots from subsequent analysis

#First, we need to tell R which spots are the highest spots on the plate
#the spots printed immediately after these can be made to have null MFIs
#This call will create a vector identifying all ref spots, plus all std spots called "Std 1"
high_targets <- c(targets_ref, grep("Std 1", row.names(annotation_targets.df)))
high_targets_disinclude <-c()

#subset high_targets to exclude high targets in the bottom of block 1
#these were printed last in both cases (reps = 1 or reps = 2)
high_targets1 <- high_targets[!between(high_targets, (index_target/2 - 12), (index_target/2))]

#To identify the spots to disinclude, we need to identify the position of the next spot in the print run.
#This is different for reps=1 or reps=2.

# For singlicate printing, the printer prints from top to bottom, right to left, with a new pickup for block 1 with different 
#targets than were just printed in block 2. 

if (reps==1){
  # If the high target is in block 2, then exclude high_target - index_target/2 
  # If the high target is in block 1, then exclude high_target + index_target/2 + 12
  high_targets_disinclude <- ifelse(high_targets1>=index_target/2, high_targets1-(index_target/2), high_targets1+(index_target/2+12))
}

# For duplicate printing, the printer prints from top to bottom, right to left, with only one pickup for both duplicate rows

if (reps==2){
  
  #subset high_targets1 to exclude high targets within (index target - 12) because these were printed last.
  high_targets2 <- high_targets1[!between(high_targets1,(index_target - 12), index_target)]
  
  #for either block 1 or block 2, exclude high target + 12
  for (i in 1:length(high_targets2)){
  high_targets_disinclude[i] <- high_targets2[i] + 12
  }
  high_targets_disinclude <- subset(high_targets_disinclude, !is.na(high_targets_disinclude))
}
remove(i)

###Convert the spots to be disincluded to NAs in background corrected data
cor.5.matrix <- cor.matrix
cor.5.matrix[high_targets_disinclude, ] <- NA

###GST subtraction!!! Do this with background corrected MFI before log transforming or anything.

#Prepare target data frame for merging (get "unique" names from annotation targets) 
target_meta2.df <- merge(target_meta.df[,1:6], annotation_targets.df, by = "Name")

#merge target data frame with data
bunny <- merge(target_meta2.df, cor.5.matrix, by.y = "row.names", by.x = "target_id_unique", sort = FALSE, all.y = TRUE)
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
plot(c(as.matrix(GST)), pch='*', col = "blue", main = "GST",
     ylab="Background Corrected MFI", xlab="GST spot", cex.main=1, cex.lab=1, cex.axis=0.7)

graphics.off()

if (reps == 1){
  GSTmean <- GST[,12:ncol(GST)]
}

if (reps == 2){
  GSTmean <- colMeans(GST[,12:ncol(GST)])
}

#subtract GST directly for each sample from all tagged targets
#set negative values to 50, which is the minimum of the background corrected data (cor.matrix)

#might have to change column number here, might need to adjust below for NAs, or remove high targets disinclude?
bundata <- bun[11:ncol(bun)]
#bundata <- bundata[,!(colnames(bundata) %in% samples_exclude1)]

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
no_tags.df <- no_tags.df[,12:ncol(no_tags.df)]

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
#Though some slides/pads may have high background, background correction may adjust this and make the data useable.
#The only way of identifying if this is true is by looking at the control spots by slide and pad.
#If the same samples have high controls, the background correction was not enough.
#Correction against plate buffers may be required...

#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The mean and cutoff will be calculated EXCLUDING any "bad" buffer targets previously set to NA.
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples

#***Samples which have deviant buffer means will be automatically excluded from further analysis
#***Buffer targets which are deviant across all pads and slides will be automatically excluded from further analysis

#Buffer assessment EXCLUDING "bad" spots: 
#Mean median corrected MFI for all data points
cor_buffer_mean <- mean(cor2.matrix[targets_buffer,], na.rm = TRUE)
#SD median corrected MFI for all data points
cor_buffer_sd <- sd(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Cut-off for deviation from mean
cor_cutoff <- cor_buffer_mean+(3*cor_buffer_sd)

###Identify deviant samples/slides/pads

# EXCLUDING "bad" spots - Generate a vector of mean buffer intensity for EACH sample (corrected person magnitude) 
cor_buffer_sample_mean <- colMeans(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Vectors to identify normal and deviant samples
cor_normal <- which(cor_buffer_sample_mean<=cor_cutoff)
cor_deviant <- which(cor_buffer_sample_mean>cor_cutoff)
#Generate a table showing which slides, pads, samples have deviant corrected buffer values
cor_sample_deviant <- samples.df[cor_deviant,]
cor_sample_deviant
write.csv(cor_sample_deviant, file = paste0(study,"_deviant_sample_buffer.csv"))

#Coefficient of variation for buffer datapoints (EXCLUDING "bad" spots), and exlcuding deviant samples
cor_buffer_cov <- cor_buffer_sd/cor_buffer_mean
cor_buffer_cov_normal <- sd(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)/mean(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)

#Plot CoV of buffer by sample - EXCLUDING "bad" spots
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

mtext(c(paste("INCLUDING Bad Buffer Spots / EXCLUDING Bad Buffer Spots:")), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Mean overall corrected buffer MFI:", round(cor_buffer_all_mean, digits=3), "/", round(cor_buffer_mean, digits=3))), side=1, cex=0.8, line=3.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("SD overall corrected buffer MFI:", round(cor_buffer_all_sd, digits=3), "/", round(cor_buffer_sd, digits=3))), side=1, cex=0.8, line=5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Cut-off overall corrected buffer MFI:", round(cor_all_cutoff, digits=3), "/", round(cor_cutoff, digits=3))), side=1, cex=0.8, line=6.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI:", round(cor_buffer_cov_all, digits=3), "/", round(cor_buffer_cov, digits=3))), side=1, cex=0.8, line=8, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI (excl. deviant samples):", round(cor_buffer_cov_all_normal, digits=3), "/", round(cor_buffer_cov_normal, digits=3))), side=1, cex=0.8, line=9.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

###Identify deviant buffer targets across all pads

#EXCLUDING "bad" spots - Mean corrected MFI for every sample for EACH target (background protein magnitude)
cor_target_mean <- rowMeans(cor2.matrix)
#Are any corrected buffer targets universally deviant, accross all pads?
deviant_buffer_targets <- Reduce(intersect, list(targets_buffer, which(cor_target_mean>cor_cutoff)))
cor_buffer_deviant.df <- annotation_targets.df[deviant_buffer_targets,]
cor_buffer_deviant.df
write.csv(cor_buffer_deviant.df, file = paste0(study,"_deviant_buffer_targets.csv"))

#list of disincluded buffer targets ("bad" spots)
removed_buffer_targets <- c()
for (i in 1:length(targets_buffer)){
  for(j in 1:length(high_targets_disinclude))
    if (targets_buffer[i] == high_targets_disinclude[j]){
      removed_buffer_targets[j] <- targets_buffer[i]
    }
}
remove(i,j)
removed_buffer_targets <- subset(removed_buffer_targets, !is.na(removed_buffer_targets))

#All slides, INCLUDING "bad" spots (ALL buffer targets)
png(filename = paste0(study, "_buffer_targets.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 3, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 1, xpd=NA, las=2)
boxplot(t(cor.matrix[targets_buffer,]),
        cex=0.5,
        ylab="Corrected MFI (log scale)", log = "y")
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(paste("Excluded buffer targets:", paste(removed_buffer_targets, collapse = ",")), las = 1)

graphics.off()

#By slide for up to 6 slides (the last slides) - INCLUDING "bad" spots (ALL buffer targets)
#mfrow is ordered by the number of rows, and columns of plots you want - so must be edited based on the number of slides
png(filename = paste0(study, "_buffer_spots_slide.tif"), width = 5, height = 4, units = "in", res = 600)
par(mfrow=c(2,3), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=2)
for (i in 1:index_slide){
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
#KG - It would be good to test this on data that has deviant buffer targets, so far they are all 0.
#KG - Actually this might be pointless as we have already calculated the values for normalization and that is the use of the buffers
cor3.matrix[deviant_buffer_targets,] <- NA

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


#### LOG TRANSFORMATION AND NORMALIZATION ###

### Log transform the data (base 2)
log.cor.matrix <- log2(cor3.matrix)

#export this matrix to compare normalized vs. not normalized data.
write.csv(t(log.cor.matrix), file = paste0(study,"_log_data.csv"))

### Normalization

###Create sample specific buffer means for normalisation
cor2_buffer_sample_mean <- colMeans(cor3.matrix[targets_buffer,], na.rm = TRUE)
cor2_buffer_sample_sd <- apply(cor3.matrix[targets_buffer,], 2, sd, na.rm = TRUE)
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
### Update this section once we have settled on names for IgG, IgA, and IgM standard curves to use in the protein key 

  #isolate data for standards, normalized and not normalized
  stds_norm <- norm.matrix[targets_std,]
  stds_pre <- log.cor.matrix[targets_std,]
  
  #check if this works to isolate the particular standard of interest - need to see names from key
  isostd_norm <- norm.matrix[c(targets_std, grep(iso, row.names(annotation_targets.df))),]

  #Plot Std 3 in Levey Jennings Style plot
  #KG - I don't know if these column numbers will hold up every time...switch to using grep?
  #average replicates for reps == 2, arithmetic mean!
  if (reps == 1){
    std_3_norm <- stds_norm[c(3),]
    std_3_pre <- stds_pre[c(3),]
    }

  if (reps == 2){
    std_3_norm <- stds_norm[c(3,9),]
    std_3_norm <- log2(apply((2^std_3_norm), 2, mean))
  
    std_3_pre <- stds_pre[c(3,9),]
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
  plot(c(std_3_norm), pch='*', col = "blue", ylim=c(min(std_3_norm, na.rm = TRUE),max(std_3_norm, na.rm=TRUE)*1.25),
     ylab="Normalized log2(MFI)", xlab="Sample (Array)")

  abline(h=std3mean)
  abline(h=std3mean+2*std3sd,lty=2)
  abline(h=std3mean-2*std3sd,lty=2)
  abline(h=std3mean+std3sd,lty=3)
  abline(h=std3mean-std3sd,lty=3) 

  plot(c(std_3_pre), pch='*', col = "darkblue", ylim=c(min(std_3_pre, na.rm = TRUE),max(std_3_pre, na.rm=TRUE)*1.25),
     ylab="log2(MFI) (NOT normalized)", xlab="Sample (Array)")

  abline(h=std3mean1)
  abline(h=std3mean1+2*std3sd1,lty=2)
  abline(h=std3mean1-2*std3sd1,lty=2)
  abline(h=std3mean1+std3sd1,lty=3)
  abline(h=std3mean1-std3sd1,lty=3)

  mtext(paste("Geometric CV, Normalized:", round(std3cv, digits=2), "%" ), side=1, cex=0.8, line=0.5, outer=TRUE, xpd=NA, adj=0)
  mtext(paste("Geometric CV, NOT Normalized:", round(std3cv1, digits=2), "%"), side=1, cex=0.8, line=1.5, outer=TRUE, xpd=NA, adj=0)

  graphics.off()

###Plot All standards together in ggplot2, but separately for rep1 and rep2
  
  #normalized
  std_norm_1 <- stds_norm[1:6,]
  std_norm_2 <- stds_norm[7:12,]
  
  std1melt <- melt(std_norm_1, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_norm_1.tif"), width = 7, height = 5, units = "in", res = 1200)
  
  ggplot(std1melt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep1 Normalized") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())
  
  graphics.off()
  
  std2melt <- melt(std_norm_2, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_norm_2.tif"), width = 7, height = 5, units = "in", res = 1200)
  
  ggplot(std2melt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep2 Normalized") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())
  
  graphics.off()
  
  #not normalized
  std_pre_1 <- stds_pre[1:6,]
  std_pre_2 <- stds_pre[7:12,]
  
  std1premelt <- melt(std_pre_1, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_pre_1.tif"), width = 7, height = 5, units = "in", res = 1200)
  
  ggplot(std1premelt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep1 Pre-Normalization") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())
  
  graphics.off()
  
  std2premelt <- melt(std_pre_2, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_pre_2.tif"), width = 7, height = 5, units = "in", res = 1200)
  
  ggplot(std2premelt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
    labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Rep2 Pre-Normalization") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())
  
  graphics.off()
  
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
  print(repR)
  
  ## Plot replicate 1 v. replicate 2 for each protein or each person and calculate correlation coefficient.
  png(filename = paste0(study, "_replicatescorrelation.tif"), width = 5, height = 4, units = "in", res = 600)
  par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
      mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
  
  plot(rep1, rep2, col="red", cex = 0.1)
  mtext(c(paste("Pearson correlation coefficient:", round(repR, digits=4))), side=3, adj=0)
  
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
for(i in 1:nrow(dup_samples)){
  sample_meta.df$exclude[which(sample_meta.df$sample_id == dup_samples$sample_id[i])] <- "yes"
}

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
  rmsamp_all <- unique(c(targets_blank, targets_buffer, targets_ref, targets_std, high_targets_disinclude))
  
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
  target_meta.df$Concentration <- as.character(target_meta.df$Concentration)
  
  #Merge with target metadata to filter based on expression tag etc.
  target.df <- merge(target_meta.df, norm_sub4.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
  
  #Save final data frame to csv file
  write.csv(norm_sub4.df, paste0(study, "_finalafterprocessing.csv"))
  
  
#Old GST subtraction code:
  # #GST subtraction - only for data with negative values set to 0 (norm4.matrix)
  # ### Subtracting Protein Tag Signal from tagged antigens - only for norm4.matrix (no negative values)
  # #Prepare data frame with GST tagged proteins only for subtraction
  # GST_antigens.df <- filter(target.df, Expression_Tag == "GST" | Expression_Tag == "GST/His")
  # GST_antigens.df <- tibble::column_to_rownames(GST_antigens.df, var="Name")
  # GST_antigens.df <- GST_antigens.df[,sapply(GST_antigens.df, is.numeric)]
  # 
  # GST <- c(grep("GST", rownames(norm_sub4.df), fixed = TRUE))
  # GST_val <- c(as.matrix(norm_sub4.df[GST,]))
  # 
  # #Plot GST values
  # png(filename = paste0(study, "_GST.tif"), width = 5, height = 3.5, units = "in", res = 1200)
  # par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(2.1,4.1,2.1,2.1))
  # plot(GST_val, pch='*', col = "blue", ylim=c(0,max(GST_val, na.rm = TRUE)*1.25), main = "GST",
  #      ylab="Normalized log2(MFI)", xlab="Sample (Array)", cex.main=1, cex.lab=1, cex.axis=0.7)
  # 
  # #print text on the plots for number of samples where GST and CD4 are above buffer
  # mtext(paste("Total Samples with GST > 0:", round(sum(norm_sub4.df[GST,] > 0, na.rm = TRUE), digits=2), 
  #             "(", round(sum(norm_sub4.df[GST,] > 0, na.rm = TRUE)/length(GST_val)*100, digits=2), "%)"), side=1, cex=0.8, line=0.5, outer=TRUE, xpd=NA, adj=0)
  # 
  # graphics.off()
  # 
  # #Subtract GST signal from GST tagged proteins
  # sub_GST_antigens.df <- data.frame(matrix(0, nrow = nrow(GST_antigens.df), ncol = ncol(GST_antigens.df)))
  # rownames(sub_GST_antigens.df) <- rownames(GST_antigens.df)
  # colnames(sub_GST_antigens.df) <- colnames(GST_antigens.df)
  # 
  # #subtract GST! There might be NA values if reps = 2.
  # for(b in 1:ncol(GST_antigens.df)){
  #   for(a in 1:nrow(GST_antigens.df)){
  #     #if GST value is NA or main value is NA, then cannot do subtraction and cannot compare with other data
  #     #so set those values to NA
  #     if(is.na(norm_sub4.df[GST,b])|is.na(GST_antigens.df[a,b])){
  #       sub_GST_antigens.df[a,b] <- NA
  #       # when the GST value is positive only, subtract GST value, 
  #       #otherwise want to leave as whatever the value was before (because GST was at or below buffer mean for that sample)
  #     } else if (norm_sub4.df[GST,b] > 0){
  #       #calculate difference in original MFI form (not log2)
  #       sub_GST_antigens.df[a,b] <- 2^GST_antigens.df[a,b] - 2^norm_sub4.df[GST,b]
  #       #can only do log2 if the difference is greater than 0, otherwise set to 0 (normalized log2 value is not above buffer)
  #       if (sub_GST_antigens.df[a,b] > 0) {
  #         sub_GST_antigens.df[a,b] <- log2(sub_GST_antigens.df[a,b])
  #         #if the log2 tag-subtracted value is negative, means normalized value is below buffer mean,
  #         #so also need to set those negatives to 0 again.
  #         if(sub_GST_antigens.df[a,b] < 0){
  #           sub_GST_antigens.df[a,b] <- 0
  #         }
  #       } else { 
  #         sub_GST_antigens.df[a,b] <- 0
  #       }
  #     } else {
  #       sub_GST_antigens.df[a,b] <- GST_antigens.df[a,b]
  #     }
  #   }
  # }
  # remove(a,b)
  # 
  # #Make another data frame where the tagged protein values are replaced by their subtracted values
  # #filter out the GST tagged targets
  # no_tags.df <- filter(target.df, !(Expression_Tag == "GST" | Expression_Tag == "GST/His"))
  # no_tags.df <- tibble::column_to_rownames(no_tags.df, var="Name")
  # no_tags.df <- no_tags.df[,sapply(no_tags.df, is.numeric)]
  # #then rbind the GST and the CD4 data frames to that one. The order of the targets
  # #shouldn't matter anymore
  # norm_sub5.df <- rbind(no_tags.df, sub_GST_antigens.df)
  # 
  # #save ALL GST subtracted data in another file
  # write.csv(norm_sub5.df, paste0(study, "_GST_subtracted_Final.csv"))
  
#Save R workspace so that can load prior to analysis 
  save.image(file= paste0(study,"_AfterProcessing.RData"))
