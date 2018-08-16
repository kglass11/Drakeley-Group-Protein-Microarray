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

### Plot standard values for each sample and assess variation - with negative normalized values

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



#Seropositivity Cutoffs

#Calculations for total number of people reactive to each Pk antigen 

#Calculations and plot for antigen breadth 


