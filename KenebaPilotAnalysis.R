#Keneba Pilot Study Data Analysis 
# Nov. 13 2018

#Post QC 

rm(list=ls())

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

workdir <- "/Users/Katie/Desktop/R files from work/Keneba pilot results/Analysis"
setwd(workdir)

#do the analysis for IgG and IgM one at a time. Load either one or the other.
#make sure to include the study name in the filenames of all plots produced.

load(file = "KenebaPi_IgGv2_AfterProcessing.RData")

###### Serum Dilution

#prepare data 
antnames <- rownames(norm_sub4.df)

tnormsub <- as.data.frame(t(norm_sub4.df))

allmeta <- merge(sample_meta_f.df,tnormsub, all.y = TRUE, by.x = "sample_id_unique", by.y = "row.names", sort = FALSE) 

#make sample_dilution a factor
allmeta$sample_dilution <- as.factor(as.character(allmeta$sample_dilution))

#make number_thaws a factor
allmeta$number_thaws <- as.factor(as.character(allmeta$number_thaws))

#separate test and controls

testmeta <- filter(allmeta, sample_type == "test" | sample_type == "test/control")
controlmeta <- filter(allmeta, sample_type == "control" | sample_type == "test/control")

#Plot all data for each antigen by each dilution as a connected line plot

#melt data frame with measure variables as antnames
testmelt <- melt(testmeta, measure.vars = antnames)
controlmelt <- melt(controlmeta, measure.vars = antnames)

#find the max value to use as the ylim of the plots so all are on same axes
max(testmelt$value, na.rm = TRUE) #8.73

#plot all data as a boxplot for each antigen all on one giant plot - test samples only
png(filename = paste0(study, "_dilution_box.tif"), width = 15, height = 6, units = "in", res = 1200)

print(ggplot(testmelt,  aes(x=variable, y = value, fill = sample_dilution)) + geom_boxplot(outlier.size = 0.2) +
        theme_bw() + labs(x = "Antigen", y = "Log2(MFI Ratio)") + 
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 8, color = "black")) +
        theme(legend.title = element_text(size = 12)) + ylim(0,9) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)))

graphics.off()

#control samples only
png(filename = paste0(study, "_dilution_box_refsera.tif"), width = 15, height = 6, units = "in", res = 1200)

print(ggplot(controlmelt,  aes(x=variable, y = value, fill = sample_dilution)) + geom_boxplot(outlier.size = 0.2) +
        theme_bw() + labs(x = "Antigen", y = "Log2(MFI Ratio)") + 
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 8, color = "black")) +
        theme(legend.title = element_text(size = 12)) + ylim(0,9) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)))

graphics.off()

#plot test samples in a line plot, separately for all antigens using facet_wrap

png(filename = paste0(study, "_dilution_line_set1.tif"), width = 13, height = 17, units = "in", res = 1200)

print(ggplot(testmelt,  aes(x=sample_dilution, y = value)) + geom_point(size = 1.2, color = "blue") +
        geom_line(aes(group = sample_id),size = 0.1) + facet_wrap(~as.factor(variable)) +
        theme_bw() + labs(x = "Dilution", y = "Log2(MFI Ratio)") + 
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 8, color = "black")) +
        theme(legend.title = element_text(size = 12)) + ylim(0,9) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) +
        theme(strip.text = element_text(face="bold", size=6), strip.background = element_rect(colour="black", size=0.3)))

graphics.off()

#repeat same plot for controls 
png(filename = paste0(study, "_dilution_line_refsera.tif"), width = 13, height = 17, units = "in", res = 1200)

print(ggplot(controlmelt,  aes(x=sample_dilution, y = value)) + geom_point(size = 1.2, color = "red") +
        geom_line(aes(group = sample_id),size = 0.1) + facet_wrap(~as.factor(variable)) +
        theme_bw() + labs(x = "Dilution", y = "Log2(MFI Ratio)") + 
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 8, color = "black")) +
        theme(legend.title = element_text(size = 12)) + ylim(0,9) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) +
        theme(strip.text = element_text(face="bold", size=6), strip.background = element_rect(colour="black", size=0.3)))

graphics.off()

tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
# ...and finally, the Paul Tol 21-color salute
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")



png(filename = paste0(study, "_dilution_line_refsera_colors.tif"), width = 13, height = 17, units = "in", res = 1200)

print(ggplot(controlmelt,  aes(x=sample_dilution, y = value, color = sample_id)) + geom_point(size = 1.2) +
        geom_line(aes(group = sample_id),size = 0.1) + facet_wrap(~as.factor(variable)) +
        theme_bw() + labs(x = "Dilution", y = "Log2(MFI Ratio)") + 
        scale_colour_manual(values = tol14rainbow) +
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 8, color = "black")) +
        theme(legend.title = element_text(size = 12)) + ylim(0,9) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) +
        theme(strip.text = element_text(face="bold", size=6), strip.background = element_rect(colour="black", size=0.3)))

graphics.off()


######## Number of freeze/thaws

levels(testmelt$number_thaws)

thaws <- filter(testmelt, number_thaws == "1" | number_thaws == "3")

#plot test samples in a line plot, separately for all antigens using facet_wrap

png(filename = paste0(study, "_n.thaws_line.tif"), width = 13, height = 17, units = "in", res = 1200)

print(ggplot(thaws,  aes(x=sample_dilution, y = value, color = sample_id, shape = number_thaws)) + geom_point(size = 1.2) +
        geom_line(aes(group = interaction(number_thaws, sample_id)),size = 0.2) + facet_wrap(~as.factor(variable)) +
        theme_bw() + labs(x = "Dilution", y = "Log2(MFI Ratio)") + 
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
        theme(axis.text = element_text(size = 10, color = "black"), legend.text = element_text(size = 8, color = "black")) +
        theme(legend.title = element_text(size = 12)) + ylim(0,9) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) +
        theme(strip.text = element_text(face="bold", size=6), strip.background = element_rect(colour="black", size=0.3)))

graphics.off()

### are we maxing out for any antigens?? 

2^max(log.cor.matrix, na.rm = TRUE)
#62603 for IgG, same as max(cor.matrix)
#44625 for IgM, same as max(cor.matrix)

#only 23 data points total for IgG data have data greater than 60,000
#out of 64512
length(which(cor.matrix > 60000))

#which antigens are they? All ref spots except one IgG high std. 
which(cor.matrix > 60000, arr.ind = TRUE)

#                             row col
# 12_REF_1                    12   9
# 504_REF_2                  504   9
# 12_REF_1                    12  10
# 252_REF_1                  252  10
# 504_REF_2                  504  10
# 504_REF_2                  504  11
# 493_REF_2                  493  25
# 264_REF_2                  264  26
# 252_REF_1                  252  35
# 12_REF_1                    12  37
# 241_REF_1                  241  37
# 469_IgG Std 1 (200ug/ml)_2 469  37
# 504_REF_2                  504  37
# 12_REF_1                    12  44
# 264_REF_2                  264  44
# 12_REF_1                    12  46
# 12_REF_1                    12  48
# 12_REF_1                    12  53
# 252_REF_1                  252  53
# 504_REF_2                  504  53
# 12_REF_1                    12  77
# 12_REF_1                    12 102
# 252_REF_1                  252 102

#### test linearity of serum dilutions?

two <- filter(allmeta, sample_dilution == "200")
four <- filter(allmeta, sample_dilution == "400")
eight <- filter(allmeta, sample_dilution == "800")

#plot the corrected factors

## calculate slope of the line from dilution 200 to 800??
  #that would be the most quantitative way to tell 
  #slope should be -1



#graphs that are antigen specific: 
# for(i in 1:length(antnames)){
#   
#   antigen = antnames[i]
#   
#   #isolate data for the antigen
#   ant1 <- filter(testmelt, variable == antigen)
#   
# }



