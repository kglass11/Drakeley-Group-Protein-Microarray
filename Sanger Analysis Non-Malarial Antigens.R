#Sanger Non-malarial Antigens Data Analysis Script
#Katie Glass
#updated: 8/23/18

#This is continuing on from the processing script in the master branch 
#called "Sanger v2 data processing.R"

#####################################
############DATA ANALYSIS############
#####################################

rm(list=ls())

#I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Data Processed
#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens"
setwd("/Users/Katie/Desktop/R files from work/100817 Sanger/Sanger Non-malarial Antigens")
getwd()

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

load(file="Sanger.2.Update.RData")
#updated to 
#older version: "SangerAfterProcessing.RData"

#I don't think we need to do this anymore
#import an updated target metadata file
# target_file <- "sanger_target_metadata_KT_v2.csv" 
# target_meta.df <- read.csv(target_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

###Prep For further analysis and FMM Cutoff determination ####

#For determining cutoffs with FMM models, using normalized data INCLUDING NEGATIVE VALUES

#Before doing any further analysis, we have to get rid of samples or targets that we are no longer 
#interested in.
#E.g. If control individuals are in our analysis, they will affect mixture model based cut-offs
#E.g. If control targets are still in our analysis, they will muck up our protein breadth estimates
#KG - for now this still includes the same protein target at different dilutions
#This means we have to subset the data, so some earlier annotations will from here on be wrong 
#(e.g. index_sample will no longer equal 96)

#Remove control protein targets
#Don't remove control samples yet, need to do tag subtraction from those samples as well, 
#and want them included in some exported data
#Do remove samples that should be excluded
norm_sub.matrix <- norm4.matrix[-rmsamp_all,(!colnames(norm4.matrix) %in% samples_exclude)]

#Replace current target names with original target names now that control targets are removed
norm_sub3.df <- merge(norm_sub.matrix, annotation_targets.df, by ="row.names", sort = FALSE)
norm_sub3.df <- tibble::column_to_rownames(norm_sub3.df, var="Row.names")
row.names(norm_sub3.df) <- norm_sub3.df$Name
norm_sub4.df <- norm_sub3.df[,1:ncol(norm_sub.matrix)]

#Make the dilution column of target_meta.df a character type
target_meta.df$Dilution <- as.character(target_meta.df$Dilution)

#Merge with target metadata to filter based on expression tag etc.
target.df <- merge(target_meta.df, norm_sub4.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)



#For this section, using normalized data with negative values set to 0 (norm4.matrix)

#Before doing any further analysis, we have to get rid of samples or targets that we are no longer 
#interested in.
#E.g. If control individuals are in our analysis, they will affect mixture model based cut-offs
#E.g. If control targets are still in our analysis, they will muck up our protein breadth estimates
#KG - for now this still includes the same protein target at different dilutions
#This means we have to subset the data, so some earlier annotations will from here on be wrong 
#(e.g. index_sample will no longer equal 96)

#Assign sample type names, to identify control and test samples (logical)
samples_test <- sample_meta.df$sample_id_unique[which(sample_meta.df$sample_type =="test")]
samples_control <- sample_meta.df$sample_id_unique[which(sample_meta.df$sample_type =="control")]

#Define a list of targets to be removed from further analysis (controls)
rmsamp_all <- unique(c(targets_blank, targets_buffer, targets_ref, targets_std, high_targets_disinclude))

#Remove control protein targets
#Don't remove control samples yet, need to do tag subtraction from those samples as well, 
#and want them included in some exported data
#Do remove samples that should be excluded
norm_sub.matrix <- norm4.matrix[-rmsamp_all,(!colnames(norm4.matrix) %in% samples_exclude)]
            
#Replace current target names with original target names now that control targets are removed
norm_sub3.df <- merge(norm_sub.matrix, annotation_targets.df, by ="row.names", sort = FALSE)
norm_sub3.df <- tibble::column_to_rownames(norm_sub3.df, var="Row.names")
row.names(norm_sub3.df) <- norm_sub3.df$Name
norm_sub4.df <- norm_sub3.df[,1:ncol(norm_sub.matrix)]

#Make the dilution column of target_meta.df a character type
target_meta.df$Dilution <- as.character(target_meta.df$Dilution)

#Merge with target metadata to filter based on expression tag etc.
target.df <- merge(target_meta.df, norm_sub4.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
  
### Subtracting Protein Tag Signal from tagged antigens - only for norm4.matrix (no negative values)
#Prepare data frame with GST tagged proteins only for subtraction
GST_antigens.df <- filter(target.df, Expression_Tag == "GST" | Expression_Tag == "GST/His")
GST_antigens.df <- tibble::column_to_rownames(GST_antigens.df, var="Name")
GST_antigens.df <- GST_antigens.df[,sapply(GST_antigens.df, is.numeric)]

#Subtract GST signal from GST tagged proteins
GST <- c(grep("GST", rownames(norm_sub4.df), fixed = TRUE))

sub_GST_antigens.df <- data.frame(matrix(0, nrow = nrow(GST_antigens.df), ncol = ncol(GST_antigens.df)))
rownames(sub_GST_antigens.df) <- rownames(GST_antigens.df)
colnames(sub_GST_antigens.df) <- colnames(GST_antigens.df)

#All the NAs have been removed already by excluding targets and samples
for(b in 1:ncol(GST_antigens.df)){
  for(a in 1:nrow(GST_antigens.df)){
  # when the GST value is positive only, subtract GST value, 
  #otherwise want to leave as whatever the value was before (because GST was at or below buffer mean for that sample)
    if (norm_sub4.df[GST,b] > 0){
        #calculate difference in original MFI form (not log2)
        sub_GST_antigens.df[a,b] <- 2^GST_antigens.df[a,b] - 2^norm_sub4.df[GST,b]
        #can only do log2 if the difference is greater than 0, otherwise set to 0 (normalized log2 value is not above buffer)
        if (sub_GST_antigens.df[a,b] > 0) {
          sub_GST_antigens.df[a,b] <- log2(sub_GST_antigens.df[a,b])
          #if the log2 tag-subtracted value is negative, means normalized value is below buffer mean,
          #so also need to set those negatives to 0 again.
            if(sub_GST_antigens.df[a,b] < 0){
              sub_GST_antigens.df[a,b] <- 0
            }
        } else { 
          sub_GST_antigens.df[a,b] <- 0
          }
    } else {
      sub_GST_antigens.df[a,b] <- GST_antigens.df[a,b]
      }
  }
}
remove(a,b)

#Prepare CD4-tagged only data frame (Sanger antigens)
CD4_antigens.df <- filter(target.df, Expression_Tag == "CD4")
CD4_antigens.df <- tibble::column_to_rownames(CD4_antigens.df, var="Name")
CD4_antigens.df <- CD4_antigens.df[,sapply(CD4_antigens.df, is.numeric)]

#Subtract CD4 from CD4 tagged antigens
CD4 <- c(grep("CD4", rownames(norm_sub4.df), fixed = TRUE))

sub_CD4_antigens.df <- data.frame(matrix(0, nrow = nrow(CD4_antigens.df), ncol = ncol(CD4_antigens.df)))
rownames(sub_CD4_antigens.df) <- rownames(CD4_antigens.df)
colnames(sub_CD4_antigens.df) <- colnames(CD4_antigens.df)

#All the NAs have been removed already by excluding targets and samples
for(b in 1:ncol(CD4_antigens.df)){
  for(a in 1:nrow(CD4_antigens.df)){
    # when the CD4 value is positive only, subtract CD4 value, 
    #otherwise want to leave as whatever the value was before (because CD4 was at or below buffer mean for that sample)
    if (norm_sub4.df[CD4,b] > 0){
      #calculate difference in original MFI form (not log2)
      sub_CD4_antigens.df[a,b] <- 2^CD4_antigens.df[a,b] - 2^norm_sub4.df[CD4,b]
      #can only do log2 if the difference is greater than 0, otherwise set to 0 (normalized log2 value is not above buffer)
      if (sub_CD4_antigens.df[a,b] > 0) {
        sub_CD4_antigens.df[a,b] <- log2(sub_CD4_antigens.df[a,b])
        #if the log2 tag-subtracted value is negative, means normalized value is below buffer mean,
        #so also need to set those negatives to 0 again.
        if(sub_CD4_antigens.df[a,b] < 0){
          sub_CD4_antigens.df[a,b] <- 0
        }
      } else { 
        sub_CD4_antigens.df[a,b] <- 0
      }
    } else {
      sub_CD4_antigens.df[a,b] <- CD4_antigens.df[a,b]
    }
  }
}
remove(a,b)

### Plot CD4 and GST signals
GST_val <- c(as.matrix(norm_sub4.df[GST,]))
CD4_val <- c(as.matrix(norm_sub4.df[CD4,]))

png(filename = paste0(study, "_Tags_GST_CD4.tif"), width = 5, height = 7.5, units = "in", res = 1200)
par(mfrow=c(2,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(GST_val, pch='*', col = "blue", ylim=c(0,max(GST_val)*1.25), main = "GST",
     ylab="Normalized log2(MFI)", xlab="Sample (Array)")

plot(CD4_val, pch='*', col = "darkblue", ylim=c(0,max(CD4_val)*1.25), main = "CD4",
     ylab="Normalized log2(MFI)", xlab="Sample (Array)")

#print text on the plots for number of samples where GST and CD4 are above buffer
mtext(paste("Total Samples with GST > 0:", round(sum(norm_sub4.df[GST,] > 0), digits=2), 
  "(", round(sum(norm_sub4.df[GST,] > 0)/length(GST_val)*100, digits=2), "%)"), side=1, cex=0.8, line=0.5, outer=TRUE, xpd=NA, adj=0)
mtext(paste("Total Samples with CD4 > 0:", round(sum(norm_sub4.df[CD4,] > 0), digits=2), 
  "(", round(sum(norm_sub4.df[CD4,] > 0)/length(CD4_val)*100, digits=2), "%)"), side=1, cex=0.8, line=1.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#Make another data frame where the tagged protein values are replaced by their subtracted values
#filter out the GST and CD4 proteins
no_tags.df <- filter(target.df, !(Expression_Tag == "CD4" | Expression_Tag == "GST" | Expression_Tag == "GST/His"))
no_tags.df <- tibble::column_to_rownames(no_tags.df, var="Name")
no_tags.df <- no_tags.df[,sapply(no_tags.df, is.numeric)]
#then rbind the GST and the CD4 data frames to that one. The order of the targets
#shouldn't matter anymore
norm_sub5.df <- rbind(no_tags.df, sub_GST_antigens.df, sub_CD4_antigens.df)

#Up to this point, everything has been the same for the malarial and non-malarial analysis.
#Save another .RData file to start here in the future. 
save.image("SangerTagSubtracted.RData")

###Seropositivity Thresholds!!!###

#DO NOT exclude any samples or control samples until AFTER seropositivity calculations!! 
#Then can subset antigens, samples, etc from the seropositivity matrix and the final data frame.

#Do seropositivity calculations on norm_sub5.df, which has excluded samples removed, but still includes controls.

###Form a seropositivity matrix based on reactivity over a cutoff derived from sample buffer background.
#The threshold is the sample buffer mean + 3SD. Then take Log2 and normalize. 
sample_cutoff <- cor2_buffer_sample_mean + 3*cor2_buffer_sample_sd
log_sample_cutoff <- log2(sample_cutoff)
norm_sample_cutoff <- log_sample_cutoff - log_buffer_sample_mean

#Tailor the norm_sample_cutoff to remove excluded samples ONLY
buffer_cutoff.matrix <- as.matrix(norm_sample_cutoff)
rownames(buffer_cutoff.matrix, colnames(norm4.matrix))
sub_cutoff <- as.matrix(buffer_cutoff.matrix[(!rownames(buffer_cutoff.matrix) %in% samples_exclude),])

#Plot the sample cutoffs for samples included in analysis
png(filename = paste0(study, "_Buffer_Cutoffs.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(sub_cutoff, pch='*', col = "blue", ylim=c(0,max(sub_cutoff)*1.25),
     ylab="Seropositivity Cutoff", xlab="Sample (Array)")

graphics.off()

#make a new cutoff where if the norm_sample_cutoff is less than 1, the cutoff is replaced with a value of 1
final_cutoff <- sub_cutoff
final_cutoff[which(final_cutoff<1)] <- 1

#How many samples have a cutoff below 1? --> 464/1361 = 34.1% (now coming out to 463/1361 = 34.01%)
sum(sub_cutoff<1)
sum(sub_cutoff<1)/length(sub_cutoff) *100

#Apply the final cutoff to ALL antigens, ALL dilutions, for test and control samples
SP_all.df <- t(apply(norm_sub5.df, 1, function(x) ((x > final_cutoff)+0)))
colnames(SP_all.df) <- colnames(norm_sub5.df)

#Make a new data frame where seropositive values will be the data and otherwise it will be NA
SP_all_data.df <- data.frame(matrix(NA, nrow = nrow(norm_sub5.df), ncol = ncol(norm_sub5.df)))
rownames(SP_all_data.df) <- rownames(norm_sub5.df)
colnames(SP_all_data.df) <- colnames(norm_sub5.df)

for(b in 1:ncol(norm_sub5.df)){
  for(a in 1:nrow(norm_sub5.df)){
    if(SP_all.df[[a,b]]==1){
      SP_all_data.df[[a,b]] <- norm_sub5.df[[a,b]] 
    }
  }
}
remove(a,b)

### Now subset both the actual data and the seropositivity information!! 

#At this point, Remove control samples for further analysis
norm_sub6.df <- norm_sub5.df[,colnames(norm_sub5.df) %in% samples_test]
SP_sub1.df <- SP_all.df[,colnames(SP_all.df) %in% samples_test]
SP_sub_data.df <- SP_all_data.df[,colnames(SP_all_data.df) %in% samples_test]

#also subset seropositivity matrices for control data
SP_con <- SP_all.df[,colnames(SP_all.df) %in% samples_control]
SP_con_data <- SP_all_data.df[,colnames(SP_all_data.df) %in% samples_control]

#Make target.df merged data frames for further use with tag-subtracted values and test samples only
target2.df <- merge(target_meta.df, norm_sub6.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
target_SP.df <- merge(target_meta.df, SP_sub1.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
target_SP_data.df <- merge(target_meta.df, SP_sub_data.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

#Also need negative control data, so make another target data frame with control samples
targetcontrol.df <- merge(target_meta.df, norm_sub5.df[,colnames(norm_sub5.df) %in% samples_control], by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
#export this to use for CP3 CV analysis
write.csv(targetcontrol.df, file = "SangerControlswithTargetinfo.csv")
#targets merged with SP control data frame
target_SP_con_data.df <- merge(target_meta.df, SP_con_data, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

#isolate non-malarial antigens
NM_target2.df <- filter(target2.df, Category == "non_malarial")
NM_targetcontrol.df <- filter(targetcontrol.df, Category == "non_malarial")

#isolate sample columns only
NMtest.df <- tibble::column_to_rownames(NM_target2.df, var="Name")
NMtest.df <- NMtest.df[,sapply(NMtest.df, is.numeric)]

NMcontrol.df <- tibble::column_to_rownames(NM_targetcontrol.df, var="Name")
NMcontrol.df <- NMcontrol.df[,sapply(NMcontrol.df, is.numeric)]

########################################
############### Plots!!! ###############
########################################

###1. Plot number of seropositive people for each NM antigen, not including negative controls

#Isolate NM antigens and test samples only for seropositivity matrix
SP_NM_test <- filter(target_SP.df, Category == "non_malarial")
SP_NM_test <- tibble::column_to_rownames(SP_NM_test, var = "Name")
SP_NM_test <- SP_NM_test[,sapply(SP_NM_test, is.numeric)]

#data frame of sums of seropositives for each antigen, sorted highest to lowest
SPpeople <- as.data.frame(as.matrix((sort(rowSums(SP_NM_test), decreasing = TRUE))))
SPpeople <- cbind(Target = rownames(SPpeople), SPpeople)
SPpeople$Target <- as.factor(SPpeople$Target)
#explicitly set factor levels to the correct order
SPpeople$Target <- factor(SPpeople$Target, levels = SPpeople$Target[order(-SPpeople$V1)])

png(filename = paste0(study, "_NM_SPpeople.tif"), width = 8, height = 4.5, units = "in", res = 1200)

ggplot(SPpeople, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + 
  ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) +
  theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black"))

graphics.off()

###2. Boxplot of seropositive data from highest to lowest median value
#Add negative controls plotted separately as a point on the same plot? or a small red line?

#isolate data for NM only (it's already test samples only)
#This data frame only has 1182 out of 1325 test samples now, 
#because the others were not seropositive to any antigens and were therefore NA, which is not numeric
SP_NM_test2 <- filter(target_SP_data.df, Category == "non_malarial")
SP_NM_test2 <- tibble::column_to_rownames(SP_NM_test2, var = "Name")
SP_NM_test2 <- SP_NM_test2[,sapply(SP_NM_test2, is.numeric)]

#melt the SP only data with Na.rm = TRUE to organize for ggplot2
SPdatamelt <- melt(as.matrix(SP_NM_test2), na.rm = TRUE)
colnames(SPdatamelt) <- c("Target", "Sample", "Normalized")

#violin plot
ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_NM_All_SP_data_violin.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_violin() + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black"))

graphics.off()

#boxplot
png(filename = paste0(study, "_NM_All_SP_data_box.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) + ylim(0,8)

graphics.off()

#add negative controls as points to the boxplot?
NM_con_data <- filter(target_SP_con_data.df, Category == "non_malarial")
NM_con_data <- tibble::column_to_rownames(NM_con_data, var = "Name")
neg_samples <-c(grep("Neg", colnames(NM_con_data)))
NM_neg_data <- NM_con_data[,neg_samples]
#There are only two values to add to the plot...one for TT and one for RubIV, both the same neg pool
#Not sure if it's worth the trouble of adding these to the plot or just mentioning in the figure caption and text.

png(filename = paste0(study, "_NM_All_SP_data_box_Neg.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) + ylim(0,8)

graphics.off()



   
#negative control data for Pf reactive targets, tag-subtracted data
neg_samples <-c(grep("Neg", colnames(target_SP_con_data.df)))
neg_data <- norm_sub5.df[-rmsamp_all, neg_samples]
Pf_neg_data <- neg_data[Pf_target_reactive==TRUE,]
Pf_neg_mean <- as.matrix(rowMeans(Pf_neg_data))













#Violin and Box Plots of data for reactive Pf antigens, sorted by highest median to lowest
ggplot(melt.Pf, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_Pf_violin_WG.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin() + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

png(filename = paste0(study, "_Pf_boxplot_WG.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

##### Repeat everything for Pf without the wheat germ antigens ######

#For seropositivity calculations, do once for Pf and once Pv, all antigens, all dilutions
Pf_no_WG.df <- filter(target2.df, Plasmodium == "Pf", !Expression_Tag=="Wheat_Germ")
Pf_no_WG.df <- tibble::column_to_rownames(Pf_no_WG.df, var="Name")
Pf_no_WG.df <- Pf_no_WG.df[,sapply(Pf_no_WG.df, is.numeric)]


#For the person_exposed calculation, only want the antigens at 1 dilution each.
#For the Sanger antigens Dilution = 0.5 
#For all other antigens, Dilution = 1
sub_Pf_no_WG.df <- filter(target2.df, (Plasmodium == "Pf" & Dilution == "1" & !Expression_Tag=="Wheat_Germ") 
                             | (Plasmodium == "Pf" & Source == "J. Rayner; WTSI" & Dilution == "0.5"))
sub_Pf_no_WG.df <- tibble::column_to_rownames(sub_Pf_no_WG.df, var="Name")
sub_Pf_no_WG.df <- sub_Pf_no_WG.df[,sapply(sub_Pf_no_WG.df, is.numeric)]

#Make seropositivity matrices for Pf , no wheat germ
SP_Pf_no_WG.df <- t(apply(Pf_no_WG.df, 1, function(x) ((x > sub_cutoff)+0)))
sub_SP_Pf_no_WG.df <- t(apply(sub_Pf_no_WG.df, 1, function(x) ((x > sub_cutoff)+0)))

###Create a threshold for overall target reactivity
#e.g. To be included in heatmaps and other analyses, perhaps targets should be reacted to by at least 5% of people?

#All Pf antigens
Pf_target_breadth2 <- rowSums(SP_Pf_no_WG.df, na.rm=TRUE)
Pf_target_reactive2 <- Pf_target_breadth2 > (ncol(SP_Pf_no_WG.df)/100)*5
cat(sum(Pf_target_reactive2), "out of", nrow(SP_Pf_no_WG.df), "Pf targets (no WG) are reactive in at least 5% of people")

###Create a threshold for overall person reactivity
#Similarly, perhaps unreactive individuals should be disincluded? Either way - this is informative.

#For person reactivity, we do not want to count multiple dilutions for each antigen 

#Sub Pf antigens
Pf_person_breadth2 <- colSums(sub_SP_Pf_no_WG.df, na.rm=TRUE)
Pf_person_exposed2 <- Pf_person_breadth2 > (nrow(sub_SP_Pf_no_WG.df)/100)*5
cat(sum(Pf_person_exposed2), "out of", ncol(sub_SP_Pf_no_WG.df), "samples are reactive to at least 5% of Pf targets")

### Export matrix of data for reactive protein targets only (cutoff mean+3SD method)

# Includes control AND test samples, no wheat germ AG
Pf.reactive.targets.matrix <- as.matrix(norm_sub5.df[Pf_target_reactive2==TRUE,])
write.csv(Pf.reactive.targets.matrix, paste0(study,"Pf_reactive_targets_data.csv")) 

###below this line, stopped making separate names for variables for w/o WG antigens.

#Pf plot of normalized data for each antigen organized by highest median
#Only include the data if the person is seropositive and an exposed person
exposed_SP_Pf.df <- SP_Pf_no_WG.df[Pf_target_reactive2==TRUE, Pf_person_exposed2==TRUE]
reactive_Pf.df <- Pf_no_WG.df[Pf_target_reactive2==TRUE,Pf_person_exposed2==TRUE]

#Export this matrix which only includes exposed individuals
write.csv(reactive_Pf.df, paste0(study,"Pf_exposed_data.csv"))

#Plot the number of seropositive people by antigen, highest to lowest for Pf
SPpeople <- as.matrix(sort(rowSums(exposed_SP_Pf.df), decreasing = TRUE))
SPpeople <- as.data.frame(SPpeople)
SPpeople <- cbind(Target = rownames(SPpeople), SPpeople)
SPpeople$Target <- as.factor(SPpeople$Target)
#explicitly set factor levels to the correct order
SPpeople$Target <- factor(SPpeople$Target, levels = SPpeople$Target[order(-SPpeople$V1)])

png(filename = paste0(study, "_Pf_num_people.tif"), width = 8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(SPpeople, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6))

graphics.off()

#Make a new data frame where seropositive values will be the number and otherwise it will be NA
SP_Pf_data.df <- data.frame(matrix(NA, nrow = nrow(reactive_Pf.df), ncol = ncol(reactive_Pf.df)))
rownames(SP_Pf_data.df) <- rownames(reactive_Pf.df)
colnames(SP_Pf_data.df) <- colnames(reactive_Pf.df)

for(b in 1:ncol(reactive_Pf.df)){
  for(a in 1:nrow(reactive_Pf.df)){
    if(exposed_SP_Pf.df[[a,b]]==1){
      SP_Pf_data.df[[a,b]] <- reactive_Pf.df[[a,b]] 
    }
  }
}
remove(a,b)

#then melt this data.frame with Na.rm = TRUE to organize for ggplot2
melt.Pf <- melt(as.matrix(SP_Pf_data.df), na.rm = TRUE)
colnames(melt.Pf) <- c("Target", "Sample", "Normalized")

#negative control data for Pf reactive targets, tag-subtracted data
neg_samples <-c(grep("Neg", colnames(norm4.matrix)))
neg_data <- norm_sub5.df[-rmsamp_all, neg_samples]
Pf_neg_data <- neg_data[Pf_target_reactive2==TRUE,]
Pf_neg_mean <- as.matrix(rowMeans(Pf_neg_data))

#Add negative control data to the plot?

#Violin and Box Plots of data for reactive Pf antigens, sorted by highest median to lowest
ggplot(melt.Pf, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_Pf_violin.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin() + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

png(filename = paste0(study, "_Pf_boxplot.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()


#######Same thing for Pv Antigens########

#Only include the data if the person is seropositive and it's an exposed individual
exposed_SP_Pv.df <- SP_Pv.df[Pv_target_reactive==TRUE, Pv_person_exposed==TRUE]
reactive_Pv.df <- Pv_antigens.df[Pv_target_reactive==TRUE,Pv_person_exposed==TRUE]

#Export this matrix which only includes exposed individuals
write.csv(reactive_Pv.df, paste0(study,"Pv_exposed_data.csv"))

#Plot the number of seropositive people by antigen, highest to lowest for Pf
SPpeoplePv <- as.matrix(sort(rowSums(exposed_SP_Pv.df), decreasing = TRUE))
SPpeoplePv <- as.data.frame(SPpeoplePv)
SPpeoplePv <- cbind(Target = rownames(SPpeoplePv), SPpeoplePv)
SPpeoplePv$Target <- as.factor(SPpeoplePv$Target)
#explicitly set factor levels to the correct order
SPpeoplePv$Target <- factor(SPpeoplePv$Target, levels = SPpeoplePv$Target[order(-SPpeoplePv$V1)])

png(filename = paste0(study, "_Pv_num_people.tif"), width = 3, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(SPpeoplePv, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6))

graphics.off()

#Make a new data frame where seropositive values will be the number and otherwise it will be NA
SP_Pv_data.df <- data.frame(matrix(NA, nrow = nrow(reactive_Pv.df), ncol = ncol(reactive_Pv.df)))
rownames(SP_Pv_data.df) <- rownames(reactive_Pv.df)
colnames(SP_Pv_data.df) <- colnames(reactive_Pv.df)

for(b in 1:ncol(reactive_Pv.df)){
  for(a in 1:nrow(reactive_Pv.df)){
    if(exposed_SP_Pv.df[[a,b]]==1){
      SP_Pv_data.df[[a,b]] <- reactive_Pv.df[[a,b]] 
    }
  }
}
remove(a,b)

#negative control data for Pv reactive targets, tag-subtracted data
Pv_neg_data <- neg_data[Pv_target_reactive==TRUE,]
Pv_neg_mean <- as.matrix(rowMeans(Pv_neg_data))

#then melt this data.frame with Na.rm = TRUE to organize for ggplot2
melt.Pv <- melt(as.matrix(SP_Pv_data.df), na.rm = TRUE)
colnames(melt.Pv) <- c("Target", "Sample", "Normalized")

#Violin and Box Plots of data for reactive Pv antigens, sorted by highest median to lowest
ggplot(melt.Pv, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_Pv_violin.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pv, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=12))

graphics.off()

png(filename = paste0(study, "_Pv_boxplot.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pv, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_boxplot(outlier.size = 0.3) + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=12))

graphics.off()


save.image("SangerAnalysis.RData")

#Old code - not using now
### Plot of geometric mean vs target, ranked from highest to lowest

# # Only using data from reactive targets and reactive test samples
# reactive.matrix <- reactive.targets.matrix[,person_exposed]
# 
# #Calculate geometric mean and geometric SD for each antigen
# #I have not done anything with the seropositivity / seronegativity here
# mean_targets <- rowMeans(reactive.matrix)
# sd_targets <- apply(reactive.matrix, 1, sd)
# 
# target_data <- data.frame(mean_targets, sd_targets, neg_mean)
# target_data <- target_data[order(-mean_targets),]
# 
# sorted_mean <- mean_targets[order(-mean_targets)]
# 
# 
# #barPlot - this still looks terrible 
# png(filename = paste0(study, "_GeoMean_Reactive_Targets.tif"), width = 5, height = 4, units = "in", res = 1200)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
# 
# barplot(sorted_mean, pch='*', col = "blue", ylim=c(0,max(mean_targets)*1.25))
# 
# axis(2, ylab="Geometric Mean Log2(Target MFI/Buffer MFI)")
# axis(1, cex.axis = 0.4, labels = as.character(row.names(target_data)), at = 1:length(sorted_mean), xlab="Antigen")
# 
# graphics.off()
# 
# ### Plot of number of seropositive individuals for each sanger antigen 
# reactive_seroposSD.matrix <- seroposSD.matrix[target_reactive==TRUE, person_exposed]
# rownames(reactive_seroposSD.matrix) <- rownames(reactive.targets.matrix)
# sub_sanger_antigens <- c(grep("(s)", rownames(reactive_seroposSD.matrix), fixed = TRUE))
# sanger_seroposSD.matrix <- reactive_seroposSD.matrix[sub_sanger_antigens,]
# 
# sanger_seroposSD.df <- tibble::rownames_to_column(sanger_seroposSD.matrix)
# 
# #Sum of people positive for each reactive sanger antigen, sorted highest to lowest
# sanger_sums <- sort(rowSums(sanger_seroposSD.matrix), decreasing = TRUE)
# #check plot in R
# barplot(sanger_sums)
# 
# #This plot still doesn't look good, I just didn't finish and make it quickly in excel to move on
# png(filename = paste0(study,"_targets_seropos_sums.tif"), width = 5.5, height = 4, units = "in", res = 600)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1), bty = "o", 
#     mgp = c(2, 0.5, 0), cex.main = 0.7, cex.axis = 0.7, cex.lab = 1, xpd=NA, las=1)
# 
# barplot(sanger_sums, xlab="Target", add=FALSE, col = "darkblue", ylim=c(0,max(sanger_sums)*1.25))
# title(main = "Number of People Seropositive for each Reactive Antigen", adj=0)
# title(ylab="Number of People", line=2.7)
# 
# graphics.off()
# 
# qplot(t(sanger_seroposSD.matrix), geom="histogram")
# 
# 
# 
# 
# 

