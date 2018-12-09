### Keneba Study Big Study Analysis 


## currently continuing on from partway through the processing script - IgG

### duplicated samples QC
#there was a replicated sample ID measured 9 times to get the interarray repeatability data 

#isolate replicated data

#post normalization and GST subtraction and replicate averaging, nothing set to 0
sample_rep <- grep("NKC821Y", colnames(normaverage.matrix))
NKC821Ydata <- normaverage.matrix[,sample_rep]

CV <- apply(NKC821Ydata, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ydata, na.rm = TRUE)
avgCV <- mean(CV, na.rm = TRUE ) * 100 #48.5% IgG

# there is one replicate where something is wrong -- > if we remove that one: 
NKC821Ydata.2 <- NKC821Ydata[,-3]

CV.2 <- apply(NKC821Ydata.2, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ydata.2, na.rm = TRUE)
avgCV.2 <- mean(CV.2, na.rm = TRUE ) * 100 #getting a negative number lol

# it doesn't make sense to get CV of negative numbers or spots where everything is close to 0.

#so let's only use the CV where the mean is > = to 1 (double buffer values)
#this is getting 60 spots total out of 144
NKC821Ygr1 <- NKC821Ydata.2[which(rowMeans(NKC821Ydata.2, na.rm = TRUE) > 1),]
CV.3 <- apply(NKC821Ygr1, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ygr1, na.rm = TRUE)
avgCV.3 <- mean(CV.3, na.rm = TRUE ) * 100 #14.3% CV IgG

# repeat without removing the deviant sample. 
NKC821Ygr2 <- NKC821Ydata[which(rowMeans(NKC821Ydata, na.rm = TRUE) > 1),]
CV.4 <- apply(NKC821Ygr2, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ygr2, na.rm = TRUE)
avgCV.4 <- mean(CV.4, na.rm = TRUE ) * 100 #27.8% CV IgG

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

