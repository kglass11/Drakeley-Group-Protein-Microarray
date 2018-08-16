#Lou Pk HCC Q-dot Study Analysis 
#August 16, 2018

#Run separately for IgG, IgA, and IgM

###Clear the environnment - OR go to Session > Clear Workspace
rm(list=ls())

#set working directory
setwd("/Users/Katie/Desktop/R files from work/Lou HCC/IgG")
getwd()

#Import data from IgG, IgA, or IgM - this script depends on importing many objects from the end of the processing scripts
load(file = "Pk_HCC_analysis_IgG.RData")

########## Additional Standard Plots ##########

###Plot All standards together in ggplot2 - for ALL THREE ISOTYPES 
Ig <- c("IgA", "Std", "IgM")

for(i in 1:length(Ig)){
  
  type = Ig[i]
  
  norm <- norm.matrix[grep(type, row.names(norm.matrix)),]
  pre <- log.cor.matrix[grep(type, row.names(log.cor.matrix)),]
  
  #normalized
  std1melt <- melt(norm, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_norm_1_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std1melt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
          labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Normalized") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
  #not normalized
  std1premelt <- melt(pre, varnames = c("Std", "Sample"))
  
  png(filename = paste0(study, "_stds_pre_1_", type, ".tif"), width = 7, height = 5, units = "in", res = 1200)
  
  print(ggplot(std1premelt, aes(x = Sample, y=value, color = Std)) + geom_point(size = 2, shape = 18) + theme_bw() +
          labs(x = "Sample", y = "Normalized Log2(MFI)", title = "Stds Pre-Normalization") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3)) +
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()))
  
  graphics.off()
  
}

###### Calculating seropositivity ######

#Seropositivity above a threshold of mean of sample specific buffer spots + 3SD. 
sample_cutoff <- cor2_buffer_sample_mean + 3*cor2_buffer_sample_sd
log_sample_cutoff <- log2(sample_cutoff)
norm_sample_cutoff <- log_sample_cutoff - log_buffer_sample_mean

#Tailor the norm_sample_cutoff to remove excluded samples and control samples
buffer_cutoff.matrix <- as.matrix(norm_sample_cutoff)
rownames(buffer_cutoff.matrix, colnames(norm4.matrix))
sub_cutoff <- buffer_cutoff.matrix[(!rownames(buffer_cutoff.matrix) %in% samples_exclude),]

#Plot the sample cutoffs for samples included in analysis
png(filename = paste0(study, "_Buffer_Cutoffs.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(sub_cutoff, pch='*', col = "blue", ylim=c(0,max(sub_cutoff)*1.25),
     ylab="Seropositivity Cutoff", xlab="Sample (Array)")

graphics.off()

#Then can apply the norm_sample_cutoff all antigens
seropos.matrix <- t(apply(norm_sub5.df, 1, function(x) ((x > sub_cutoff)+0)))

#Select Pk PCR+ samples from SEROPOSITIVITY MATRIX
seroposT <- t(seropos.matrix)

SP.meta <- merge(sample_meta_f.df, seroposT, by.y = "row.names", by.x = "sample_id_unique", sort = FALSE)
SP.Pk <- filter(SP.meta, pcr == "Pk")
rownames(SP.Pk) <- SP.Pk$sample_id_unique

SP.Pk.data <- SP.Pk[,(ncol(sample_meta_f.df)+1):ncol(SP.Pk)]

#Select Lou antigens 
SP.Pk.meta <- merge(target_meta.df, t(SP.Pk.data), by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
SP.Pk.Lou.meta <- filter(onlySPmeta.df, Source == "Lou")
SP.Pk.Lou <- tibble::column_to_rownames(SP.Pk.Lou, var="Name")
SP.Pk.Lou <- SP.Pk.Lou[,ncol(target_meta.df):ncol(SP.Pk.Lou)]

#Make a data frame with only seropositive data, NA for everything else
#This is for ALL SAMPLES and ANTIGENS in case you want the other data later
onlySP.df <- as.data.frame(matrix(NA, nrow = nrow(norm_sub5.df), ncol = ncol(norm_sub5.df)))
rownames(onlySP.df) <- rownames(norm_sub5.df)
colnames(onlySP.df) <- colnames(norm_sub5.df)

for(i in 1:nrow(norm_sub5.df)){
  for(k in 1:ncol(norm_sub5.df)){
    if (seropos.matrix[i,k] == 1){
      onlySP.df[i,k] <- norm_sub5.df[i,k]
    }
  }
}

#seropositive data only, filtering out Lou antigens
onlySPmeta.df <- merge(target_meta.df, onlySP.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

onlySPmetaLou.df <- filter(onlySPmeta.df, Source == "Lou")
tacos.df <- tibble::column_to_rownames(onlySPmetaLou.df, var="Name")
tacos.df <- tacos.df[,ncol(target_meta.df):ncol(tacos.df)]

#now select only Pk samples
tacos.meta <- merge(sample_meta_f.df, t(tacos.df), by.y = "row.names", by.x = "sample_id_unique", sort = FALSE)
tacos.Pk <- filter(tacos.meta, pcr == "Pk")
rownames(tacos.Pk) <- tacos.Pk$sample_id_unique

tacos.Pk.data <- tacos.Pk[,(ncol(sample_meta_f.df)+1):ncol(tacos.Pk)]

############ Calculations for total number of people reactive to each Pk antigen ##########

############ Calculations and plot for antigen breadth #############




