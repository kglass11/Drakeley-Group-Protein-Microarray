#Sanger Data Analysis Script
#Katie Glass
#updated: 9/12/18

#####################################
############DATA ANALYSIS############
#####################################

rm(list=ls())

#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens/Sanger NM V2"
#"I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Non-malarial Antigens"
#
setwd("/Users/Katie/Desktop/R files from work/100817 Sanger/Sanger NM v3")
getwd()

require("gtools")

library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(reshape2)
library(corrplot)
library(ggbeeswarm)

load(file="Sanger.2.Update.RData")
load(file = "sangerNMcutoffsfinal.RData")
#note - the sangerNMcutoffsfinal.RData is not there, cannot open...find this later.

###### Isolate data for non-malarial antigens and test samples only (no controls) or test samples and controls
    ###negative values are included here
    ###this first section is repeated from FMM script (Sanger Analysis Non-Malarial Antigens)

#Remove control protein targets
#Remove samples that should be excluded and isolate test samples - leaves 1325 samples
#This means that the two neg pools and 6 malaria positive controls are not affecting the FMMs
subdata <- norm3.matrix[-rmsamp_all,!(colnames(norm3.matrix) %in% samples_exclude)]
testdata <- subdata[,(colnames(subdata) %in% samples_test)]

#Replace current target names with original target names now that control targets are removed
testdata.df <- merge(testdata, annotation_targets.df, by ="row.names", sort = FALSE)
testdata.df <- tibble::column_to_rownames(testdata.df, var="Row.names")
row.names(testdata.df) <- testdata.df$Name
Ftestdata <- testdata.df[,1:ncol(testdata)]

#do the same thing for the matrix that has the controls - need this later
#Replace current target names with original target names now that control targets are removed
subdata.df <- merge(subdata, annotation_targets.df, by ="row.names", sort = FALSE)
subdata.df <- tibble::column_to_rownames(subdata.df, var="Row.names")
row.names(subdata.df) <- subdata.df$Name
Fsubdata <- subdata.df[,1:ncol(subdata)]

#Merge with target metadata to filter based on expression tag etc.
target.df <- merge(target_meta.df, Ftestdata, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
targetsub.df <- merge(target_meta.df, Fsubdata, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

#isolate data for non-malarial antigens
NMdatameta <- filter(target.df, Category == "non_malarial")
NMdata <- NMdatameta[,(ncol(target_meta.df)+1):ncol(NMdatameta)]
rownames(NMdata) <- NMdatameta$Name

min(NMdata) #-7.300968
max(NMdata) #8.524411

#repeat for subdata (test and control samples)
NMallsampmeta <- filter(targetsub.df, Category == "non_malarial")
NMallsampdata <- NMallsampmeta[,(ncol(target_meta.df)+1):ncol(NMallsampmeta)]
rownames(NMallsampdata) <- NMallsampmeta$Name

#Transpose data so antigens are columns - it's now a matrix
tNMdata <- t(NMdata)

tNMallsampdata <- t(NMallsampdata)

##### Remove antigens which we are not analyzing from cutoff matrix and data matrix
rmant <- c("CT706","Pertussis JNIH-5 [100] *", "CT110", "Pertussis JNIH-3 [1] *", "Pertussis JNIH-3 [0.1] *")

finalcut <- as.matrix(cutoffsavedfinal[!rownames(cutoffsavedfinal) %in% rmant,])

itta <- as.data.frame(tNMdata[,!colnames(tNMdata) %in% rmant])

############### Correlation between antigens #################

### do correlation plots with ALL data, including negative values. 

#run correlation significance test with pearson's method
restime <- cor.mtest(itta, conf.level = .95, na.rm = TRUE)

#Plot with significance showing - x over NOT significant
png(filename = paste0("Sig_allANT_Correlogram.tif"), width = 10, height = 9.5, units = "in", res = 1200)

print(corrplot.mixed(cor(itta, use = "complete.obs"), p.mat = as.matrix(restime$p), sig.level = .05, tl.col="black", order = "FPC", 
    tl.pos = "lt", tl.cex = 0.5, number.cex = 0.7))

graphics.off()

#Plot without significance showing - only a few not significant because of high n. the significance test is basically pointless. 
#note, the above plot is in a different order than this one because of the significance not working with alphabetical order
png(filename = paste0("allANT_Correlogram.tif"), width = 10, height = 9.5, units = "in", res = 1200)

print(corrplot.mixed(cor(itta, use = "complete.obs"), tl.col="black", order = "alphabet", 
                     tl.pos = "lt", tl.cex = 0.7, number.cex = 0.75))

graphics.off()

#Scatter plot for HSV1 vs HSV2
png(filename = paste0(study, "_HSV1_v_HSV2.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

plot(itta$HSV1, itta$HSV2, col="red", cex = 0.1, xlab = "HSV1", ylab = "HSV2")
abline(h=1.3977643) #cutoff for HSV2
abline(v=1.4123652) #cutoff for HSV1

graphics.off()

################### Clustering of all samples with all antigens ###########################

### All data, including negative values
# ? Include a couple malarial antigens as well? Haven't yet, but good idea

### Add dengue PC1 instead of dengue 1-4

ittacluster1 <- cbind(denguedata2[,5], itta)

#column 1 = dengue PCA PC1

### remove repeated antigens and dengue 1-4
colnames(ittacluster1)

rmant2 <- c("Tg", "DENV1-NS1","DENV2-NS1","Pertussis JNIH-5 [0.1] *","Pertussis JNIH-5 [10] *",
            "DENV3-NS1","DENV4-NS1","Pertussis JNIH-3 [10] *","Pertussis JNIH-5 [1] *")
ittacluster1 <- ittacluster1[,!colnames(ittacluster1) %in% rmant2]

#scale the data because principle component not on same scale 
#this definitely changed the cluster analysis
ittacluster <- as.data.frame(scale(ittacluster1))

#### Hierarchical clustering
#need to supply a distance matrix, which is the distance of every point to every other point
d <- dist(ittacluster)

#usually need to try different algorithms, ward.D2 pre-selected dunno why though
fitH <- hclust(d, "ward.D2")
plot(fitH)
rect.hclust(fitH, k = 5, border = "red")

hclusters <- cutree(fitH, k = 5)
hclusters

#need to see if the clusters mean anything --> do they group with age etc

#save this as a large plot and see if any antigens show the clusters well
png(filename = paste0(study, "HclustPlotScatter.tif"), width = 17, height = 13, units = "in", res = 600)
plot(ittacluster, col = hclusters)

graphics.off()
#nothing really jumping out clusters are all super overlapping - try a heatmap

ittaclusterH <- as.data.frame(cbind(hclusters, ittacluster))

ittaclusterH$hclusters <- as.factor(ittaclusterH$hclusters)

#plot with stat ellipse which shows 95% confidence interval
ggplot(ittaclusterH, aes(x=ittaclusterH[,2], y = ittaclusterH$PgP3, color = ittaclusterH$hclusters, fill = ittaclusterH$hclusters)) + 
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, color = "black")

#heatmap - didn't work as is, need to get the dendrogram object out
heatmap.2(as.matrix(ittacluster), scale = "none", )



#do other clustering methods, then add all the clusters to ittaclusterH

#then merge with metadata to see if clusters relate to anything



length(ittaclusterH$hclusters)












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

########################################
############### Plots!!! ###############
########################################

###1. Plot number of seropositive people for each NM antigen, not including negative controls
#there are 1325 test samples
#there are 1361 with all controls. There are only two negative control pools. 

SP_test.df <- SP_all.df[(rownames(SP_all.df) %in% samples_test),]
dim(SP_test.df)

#data frame of sums of seropositives for each antigen, sorted highest to lowest
SPpeople <- as.data.frame(as.matrix((sort(rowSums(t(SP_test.df)), decreasing = TRUE))))
SPpeople <- cbind(Target = rownames(SPpeople), SPpeople)
SPpeople$Target <- as.factor(SPpeople$Target)
#explicitly set factor levels to the correct order
SPpeople$Target <- factor(SPpeople$Target, levels = SPpeople$Target[order(-SPpeople$V1)])

png(filename = paste0(study, "_NM_SPpeople.tif"), width = 8, height = 4.5, units = "in", res = 1200)

ggplot(SPpeople, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + 
  ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8)) +
  theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black"))

graphics.off()

###2. Boxplot of seropositive data from highest to lowest median value
#Add negative controls plotted separately as a point on the same plot? or a small red line?

SP_NM_test2 <- SP_all_data.df[(rownames(SP_all_data.df) %in% samples_test),]
dim(SP_NM_test2)

#melt the SP only data with Na.rm = TRUE to organize for ggplot2
SPdatamelt <- melt(as.matrix(SP_NM_test2), na.rm = TRUE)
colnames(SPdatamelt) <- c("Sample","Target", "Normalized")

#violin plot
ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_NM_All_SP_data_violin.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_violin(fill="lightblue", scale = "width") + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black"))

graphics.off()

#boxplot vertical
png(filename = paste0(study, "_NM_All_SP_data_box.tif"), width = 5, height = 8, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_boxplot(outlier.size = 0.3, fill ="lightblue") + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) + ylim(0,8)

graphics.off()

#boxplot horizontal
png(filename = paste0(study, "_NM_All_SP_data_box_H.tif"), width = 8, height = 5, units = "in", res = 1200)

ggplot(SPdatamelt, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + theme_bw() +
  geom_boxplot(outlier.size = 0.3, fill ="lightblue") + xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10)) + theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
  theme(axis.text.x = element_text(color = "black")) + ylim(0,8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8))

graphics.off()

#add negative controls as points to the boxplot?
neg_samples <-c(grep("Neg", rownames(SP_all_data.df)))
NM_neg_data <- SP_all_data.df[neg_samples,]

neg_data_melt <- melt(NM_neg_data, na.rm = TRUE)
neg_data_melt

#     variable  value
#6       CMV 2.305173
#45       TT 3.542611
#46       TT 3.716070

#There are only three values to add to the plot.
#Not worth adding these to the plot - just mention them in the figure caption and text.

###3. Plot Antigen Continuous data by age, gender, occupation, hospitalized by malaria, 
  #ethnicity, number of episodes of malaria in life

sample_meta_f.df$AgeBin <- NA

#add age bins to sample_meta_f.df 
for(i in 1:nrow(sample_meta_f.df)){
  if (is.na(sample_meta_f.df$"Age in years"[i])) {i = i +1}
  else if (sample_meta_f.df$"Age in years"[i] < 5) {sample_meta_f.df$AgeBin[i] <- "< 5"}
  else if (sample_meta_f.df$"Age in years"[i] >= 5 & sample_meta_f.df$"Age in years"[i] < 15) {sample_meta_f.df$AgeBin[i] <- "5-14"}
  else if (sample_meta_f.df$"Age in years"[i] >= 15 & sample_meta_f.df$"Age in years"[i] < 25) {sample_meta_f.df$AgeBin[i] <- "15-24"}
  else if (sample_meta_f.df$"Age in years"[i] >= 25 & sample_meta_f.df$"Age in years"[i] < 35) {sample_meta_f.df$AgeBin[i] <- "25-34"}
  else if (sample_meta_f.df$"Age in years"[i] >= 35 & sample_meta_f.df$"Age in years"[i] < 50) {sample_meta_f.df$AgeBin[i] <- "35-49"}
  else if (sample_meta_f.df$"Age in years"[i] >= 50 & sample_meta_f.df$"Age in years"[i] < 70) {sample_meta_f.df$AgeBin[i] <- "50-69"}
  else if (sample_meta_f.df$"Age in years"[i] >= 70 & sample_meta_f.df$"Age in years"[i] < 100) {sample_meta_f.df$AgeBin[i] <- "70-100"}
}

#explicitly set factor levels for age bins
sample_meta_f.df$AgeBin <- factor(sample_meta_f.df$AgeBin, levels = c("< 5", "5-14","15-24","25-34", "35-49", "50-69", "70-100"))


#merge relevant data with sample meta data 
subtacos <- merge(sample_meta_f.df, SP_NM_test2 , by.x = "sample_id_unique", by.y = "row.names", sort = FALSE)

antnames <- colnames(SP_NM_test2)

#subtacos$day <- factor(subtacos$day, levels = as.character(c("0", "7", "28")))

#need to melt the data for ggplot2

subtacosm <- melt(subtacos, measure.vars = antnames, na.rm = TRUE)

subtacosm$"Age in years" <- as.numeric(subtacosm$"Age in years")
subtacosm$"Number of episodes of malaria in life" <- as.numeric(subtacosm$"Number of episodes of malaria in life")

#plot with selected epi variables vs antibody response

for(i in 1:length(antnames)){
  
  antigen = antnames[i]
  
  #isolate data for the antigen
  ant1 <- filter(subtacosm, variable == antigen)
  
  #isolate cutoff for the antigen
  cut1 <- finalcut[i]
  
  #age in years (no bins) scatter plot
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.age.tif"), width = 3.5, height = 3, units = "in", res = 1200)
  
  print(ggplot(ant1, aes(x = ant1$"Age in years", y = value)) + geom_point(color = "blue", shape = 17, size = 0.5) +
          theme_bw() + labs(x = "Age", y = "Log2(MFI Ratio)", title = antigen) + 
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12))+ xlim(0,100) + ylim(0,8) +
          geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  
  graphics.off()
  
  #age bins bee swarm and violin plot -- change the plot below! 
  ant1bin <- filter(ant1, !AgeBin == "NA")
  
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.ageBINS_V_bee.tif"), width = 7, height = 4, units = "in", res = 1200)
  
  print(ggplot(ant1bin, aes(x = AgeBin, y = value, color = AgeBin)) + geom_violin(scale = "width", color = "black") +
          theme_bw() + labs(x = "Age", y = "Log2(MFI Ratio)", title = antigen) + geom_beeswarm(cex = 0.75, size = 0.5, show.legend = F) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12)) + ylim(0,8) +
          geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  
  graphics.off()
  
  # by gender boxplot
  ant1gender <- filter(ant1, Gender == "male" | Gender == "female")
  
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.Gender.tif"), width = 2.7, height = 3, units = "in", res = 1200)
  
  print(ggplot(ant1gender, aes(x = Gender, y = value, fill = Gender)) + geom_boxplot(outlier.size = 0.3, show.legend=F) +
          theme_bw() + labs(x = "Gender", y = "Log2(MFI Ratio)", title = antigen) + 
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12)) + ylim(0,8) +
          geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  
  graphics.off()
  
  #by gender beeswarm and violin plot
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.Gender_V_bee.tif"), width = 3, height = 4, units = "in", res = 1200)
  
  print(ggplot(ant1gender, aes(x = Gender, y = value, color = Gender)) + geom_violin(scale = "width", color = "black") +
          theme_bw() + labs(x = "Gender", y = "Log2(MFI Ratio)", title = antigen) + geom_beeswarm(cex = 1, size = 0.5, show.legend = F) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12)) + ylim(0,8) +
          geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  
  graphics.off()
  
  #by occupation beeswarm and violin plot
  ant1occ <- filter(ant1, !(Occupation == "" | Occupation == "NA"))
  
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.Occup_V_bee.tif"), width = 7, height = 4, units = "in", res = 1200)
  
  print(ggplot(ant1occ, aes(x = Occupation, y = value, color = Occupation)) + geom_violin(scale = "width", color = "black") +
          theme_bw() + labs(x = "Occupation", y = "Log2(MFI Ratio)", title = antigen) + geom_beeswarm(cex = .3, size = 0.5, show.legend = F) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12)) + ylim(0,8) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 9)) +
          geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  
  graphics.off()
  
  #by hospitalized by malaria beeswarm and violin plot 
  ant1hosp <- filter(ant1, ant1$"Hospitalised by malaria" == "yes" |  ant1$"Hospitalised by malaria" == "no")
  
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.HospM_V_bee.tif"), width = 3, height = 4, units = "in", res = 1200)
  
  print(ggplot(ant1hosp, aes(x = ant1hosp$"Hospitalised by malaria", y = value, color = ant1hosp$"Hospitalised by malaria")) + geom_violin(scale = "width", color = "black") +
          theme_bw() + labs(x = "Hospitalized by Malaria", y = "Log2(MFI Ratio)", title = antigen) + geom_beeswarm(cex = 1, size = 0.5, show.legend = F) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12)) + ylim(0,8) +
          geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  
  graphics.off()
  
  #by number malaria episodes scatter plot - this is probably not worth doing, but can tell if it is from this plot
  ant1malep <- filter(ant1, !(ant1$"Number of episodes of malaria in life" == "" | ant1$"Number of episodes of malaria in life" == "NA"))
  
  png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.MalEp.tif"), width = 3.5, height = 3, units = "in", res = 1200)
  
  print(ggplot(ant1malep, aes(x = ant1malep$"Number of episodes of malaria in life", y = value)) + geom_point(color = "blue", shape = 17, size = 0.5) +
          theme_bw() + labs(x = "Number of malaria episodes", y = "Log2(MFI Ratio)", title = antigen) + 
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
          theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
          theme(legend.title = element_text(size = 12))+ xlim(0,10) + ylim(0,8) +
          geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  
  graphics.off()
  
  #by number malaria episodes beeswarm and violin plot -- this looks stupid, violin part not working  
  # 
  # 
  # png(filename = paste0(study, "_", antigen,"_SP_Ab_vs.MalEp_V_bee.tif"), width = 7, height = 4, units = "in", res = 1200)
  # 
  # print(ggplot(ant1malep, aes(x = ant1malep$"Number of episodes of malaria in life", y = value, color = ant1malep$"Number of episodes of malaria in life")) + geom_violin(scale = "width", color = "black") +
  #         theme_bw() + labs(x = "Number of malaria episodes", y = "Log2(MFI Ratio)", title = antigen) + geom_beeswarm(cex = .3, size = 0.5, show.legend = F) +
  #         theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
  #         theme(axis.text = element_text(size = 12, color = "black"), legend.text = element_text(size = 12, color = "black")) +
  #         theme(legend.title = element_text(size = 12)) + ylim(0,8) +
  #         geom_hline(yintercept=cut1, linetype="dashed", color = "black", size=0.2))
  # 
  # graphics.off()
  
}
