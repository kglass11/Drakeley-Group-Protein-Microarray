### Keneba Study Big Study Analysis 


## currently continuing on from partway through the processing script - IgG

### duplicated samples QC
#there was a replicated sample ID measured 9 times to get the interarray repeatability data 

#isolate replicated data

#post normalization and GST subtraction and replicate averaging, nothing set to 0
sample_rep <- grep("NKC821Y", colnames(normaverage.matrix))
NKC821Ydata <- normaverage.matrix[,sample_rep]

CV <- apply(NKC821Ydata, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ydata, na.rm = TRUE)
avgCV <- mean(CV, na.rm = TRUE ) * 100 #48.5% IgG, 40.7 IgM

# there is one replicate where something is wrong -- > if we remove that one: 
NKC821Ydata.2 <- NKC821Ydata[,-3]

CV.2 <- apply(NKC821Ydata.2, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ydata.2, na.rm = TRUE)
avgCV.2 <- mean(CV.2, na.rm = TRUE ) * 100 #getting a negative number lol, 8.88% for IgM

# it doesn't make sense to get CV of negative numbers or spots where everything is close to 0.

#so let's only use the CV where the mean is > = to 1 (double buffer values)
#this is getting 60 spots total out of 144
NKC821Ygr1 <- NKC821Ydata.2[which(rowMeans(NKC821Ydata.2, na.rm = TRUE) > 1),]
CV.3 <- apply(NKC821Ygr1, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ygr1, na.rm = TRUE)
avgCV.3 <- mean(CV.3, na.rm = TRUE ) * 100 #14.3% CV IgG, 11.7% for IgM

# repeat without removing the deviant sample. 
NKC821Ygr2 <- NKC821Ydata[which(rowMeans(NKC821Ydata, na.rm = TRUE) > 1),]
CV.4 <- apply(NKC821Ygr2, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ygr2, na.rm = TRUE)
avgCV.4 <- mean(CV.4, na.rm = TRUE ) * 100 #27.8% CV IgG, 16.5% for IgM

#trying to see what's up with that sample -- look at the replicates
#everything fine with those - all 4 replicates are bad for that sample and correlate very closely
NKC821YdataALL <- norm.matrix[,sample_rep]

#look at buffer means for normalization - that one is only slightly higher
#range 6.58 - 6.95, doesn't explain the massive difference in signal 
log_buffer_sample_mean[grep("NKC821Y", names(log_buffer_sample_mean))]

# check if there was weirdly high background / if the problem was there in the raw data
#the problem is present in the original data (fore) and there is no high background
NKC821Yback <- back.matrix[,grep("NKC821Y", colnames(back.matrix))]
NKC821Yfor <- fore.matrix[,grep("NKC821Y", colnames(fore.matrix))]

####### IgG vs IgM standard reactivity 

IgG_high_mean <- mean(normaverage.matrix[grep("IgG Std 1", rownames(normaverage.matrix)),], na.rm = TRUE)
IgM_high_mean <- mean(normaverage.matrix[grep("IgM Std 1", rownames(normaverage.matrix)),], na.rm = TRUE)
IgMfc_high_mean <- mean(normaverage.matrix[grep("IgM Fc Std 1", rownames(normaverage.matrix)),], na.rm = TRUE)

if(iso == "IgG"){
  stdratio <- IgG_high_mean/IgM_high_mean
  eratio <- 2^stdratio
  stdratioFc <- IgG_high_mean/IgMfc_high_mean
  eratioFc <- 2^stdratioFc
}

if(iso == "IgM"){
  stdratio <- IgM_high_mean/IgG_high_mean
  eratio <- 2^stdratio
  stdratioFc <- IgMfc_high_mean/IgG_high_mean
  eratioFc <- 2^stdratioFc
}

