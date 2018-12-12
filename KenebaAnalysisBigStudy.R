### Keneba Study Big Study Analysis 


## currently continuing on from partway through the processing script - IgG

### duplicated samples QC
#there was a replicated sample ID measured 9 times to get the interarray repeatability data 

#isolate replicated data

#post normalization and GST subtraction and replicate averaging, nothing set to 0
sample_rep <- grep("NKC821Y", colnames(quadrulemat))
NKC821Ydata <- quadrulemat[,sample_rep]

CV <- apply(NKC821Ydata, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ydata, na.rm = TRUE)
avgCV <- mean(CV, na.rm = TRUE ) * 100 

# but it doesn't make sense to get CV of negative numbers or spots where everything is close to 0

#so let's only use the CV where the mean is > = to 1 (double buffer values)
#this is getting 60 spots total out of 144
NKC821Ygr1 <- NKC821Ydata[which(rowMeans(NKC821Ydata, na.rm = TRUE) > 1),]
CV.3 <- apply(NKC821Ygr1, 1, sd, na.rm = TRUE) / rowMeans(NKC821Ygr1, na.rm = TRUE)
avgCV.3 <- mean(CV.3, na.rm = TRUE ) * 100 #12.15% CV IgG, 11.7% for IgM

####### IgG vs IgM standard reactivity 

IgG_high_mean <- mean(quadrulemat[grep("IgG Std 1", rownames(quadrulemat)),], na.rm = TRUE)
IgM_high_mean <- mean(quadrulemat[grep("IgM Std 1", rownames(quadrulemat)),], na.rm = TRUE)
IgMfc_high_mean <- mean(quadrulemat[grep("IgM Fc Std 1", rownames(quadrulemat)),], na.rm = TRUE)

if(iso == "IgG"){
  stdratio <- IgG_high_mean/IgM_high_mean
  eratio <- 2^stdratio #7.945517
  stdratioFc <- IgG_high_mean/IgMfc_high_mean
  eratioFc <- 2^stdratioFc #4.89357e+40
}

if(iso == "IgM"){
  stdratio <- IgM_high_mean/IgG_high_mean
  eratio <- 2^stdratio
  stdratioFc <- IgMfc_high_mean/IgG_high_mean
  eratioFc <- 2^stdratioFc
}

