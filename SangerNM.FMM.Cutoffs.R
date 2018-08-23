#Sanger NM FMM Seropositivity Cutoffs

rm(list=ls())

setwd("/Users/Katie/Desktop/R files from work/100817 Sanger/Sanger Data Processing KG")

load(file = "Sanger.2.Update.RData")

require("gtools")

library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)

#define functions that take results of mixture model to set cutoffs

f<-function(x,prob,lambda,mu,sigma,k,k1){
  
  lista<-order(mu)
  
  lambda<-lambda[lista]
  mu<-mu[lista]	
  sigma<-sigma[lista]	
  
  p<-vector('numeric',k)
  
  for(i in 1:k)p[i]<-lambda[i]*dnorm(x,mu[i],sigma[i])
  
  prob-sum(p[1:k1])/sum(p)
  
}

f2<-function(x,prob,lambda,mu,sigma,k,k1){
  
  lista<-order(mu)
  
  lambda<-lambda[lista]
  mu<-mu[lista]	
  sigma<-sigma[lista]	
  
  p<-vector('numeric',k)
  
  for(i in 1:k)p[i]<-lambda[i]*dnorm(x,mu[i],sigma[i])
  
  prob-sum(p[k1:k])/sum(p)
  
}

#2 Separate data to antigen of interest

#Yuyen's selection of his data:
#pos.data3<-pos.data2[is.na(pos.data2$"wb12376")==F&pos.data2$"wb12376">0&pos.data2$age>=1,]

#For this study we want to start with GST subtracted data, 
#without negative normalized values set to 0 --> start with trans.norm.df

#remove control samples and only use test samples in the future?
antibody <- as.numeric(c(trans.norm.df$"53_Tg_1"))

min(antibody,na.rm=TRUE)

max(antibody,na.rm=TRUE)

mean(antibody,na.rm=TRUE)

sd(antibody,na.rm=TRUE)

sum(is.na(antibody)==T)

antibody1<-sort(antibody)

#4 Plot Mixture Model

png(filename = paste0(study,"_Pre_FMM_53_Tg.tif"), width = 8, height = 4, units = "in", res = 600)
par(mfrow=c(1,2))

plot(density(antibody1),xlab='Log2(MFI Ratio)',main='')
title('A',adj=0,cex.main=1.5)
title('Tg',adj=0.5)

qqnorm(antibody1,las=1,pch=21,bg='grey',cex=0.75)
qqline(antibody1)
title('B',adj=0,cex.main=1.5)

graphics.off()

#Remember to change mu based on the density plot
fit.ab2<-normalmixEM(antibody1,lambda=c(0.5,0.5),mu=c(0,5),k=2)

#5 Plot Cut off Value
summary(fit.ab2)

#cutoff below which is negative - manually set interval end points
cutoff<-uniroot(f,c(-3,9),prob=0.90,lambda=fit.ab2$lambda,mu=fit.ab2$mu,sigma=fit.ab2$sigma,k=2,k1=1)$root

#cutoff above which is positive - manually set interval end points
cutoff2<-uniroot(f2,c(-3,9),prob=0.90, lambda=fit.ab2$lambda,mu=fit.ab2$mu,sigma=fit.ab2$sigma,k=2,k1=2)$root

#percentage positive 
pos <- sum(antibody1>cutoff2)/length(antibody1) * 100

#percentage negative
neg <- sum(antibody1<cutoff)/length(antibody1) * 100

#percentage indeterminate
indet <- (1-sum(antibody1>cutoff2)/length(antibody1)-sum(antibody1<cutoff)/length(antibody1)) *100

#Plot cutoffs
png(filename = paste0(study,"_FMM_Cutoffs_53_Tg.tif"), width = 5, height = 5, units = "in", res = 600)
par(mfrow = c(1,1), mar = c(5, 5, 2, 2), oma = c(6, 1, 1, 1))

plot(antibody1,fit.ab2$posterior[,1],type='n',xlim=c(0,10),lwd=2,ylim=c(0,1),col='green',las=1,xlab='Log2(MFI Ratio)',ylab='classification probability')

rect(cutoff,-0.04,cutoff2,1.04,col='light grey',lwd=1.5)
title('C',adj=0,cex.main=1.5)
lines(antibody1,fit.ab2$posterior[,2],lwd=2,col='red')

abline(h=0.90,lwd=1.5,lty=2)
abline(v=cutoff,lwd=1)
abline(v=cutoff2,lwd=1)

lines(antibody1,fit.ab2$posterior[,1],lwd=2,col='blue')
lines(antibody1,fit.ab2$posterior[,2],lwd=2,col='purple')
legend(1125,0.8,c(expression(S^'-'),expression(S^'+')),lty=c(1,1),lwd=2,col=c('blue','purple'))

mtext(c(paste("Positive Cutoff:", round(cutoff2, digits=3), "; ", round(pos, digits = 2), "% of samples")), side=1, cex=0.8, line=1.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Negative Cutoff:", round(cutoff, digits=3),"; ", round(neg, digits = 2), "% of samples")), side=1, cex=0.8, line=3, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste(round(indet, digits=2),"% of samples are indeterminate")), side=1, cex=0.8, line=4.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#Lindsey's code / plots (adapted from her function "cutoff": 

plot.data <- antibody1

fit <- normalmixEM(na.omit(plot.data),2)
print("FMM 2 component (MFI)")
summary(fit)

mu1 <- fit$mu[1]
mu2 <- fit$mu[2]
sig1 <- fit$sigma[1]
sig2 <- fit$sigma[2]

min_comp1 <- which(fit$mu == min(fit$mu))

cut <- cut_all[i] <- fit$mu[min_comp1]+no_sd*sqrt(fit$sigma[min_comp1])

plot_x <- seq(0,max(plot.data,na.rm=T),1)

gauss1a <- dnorm(plot_x,mu1,sig1)
gauss2a <- dnorm(plot_x,mu2,sig2)

png(filename = paste0(study,"_FMM_Cutoffs_53_Tg_Lindsey.tif"), width = 5, height = 5, units = "in", res = 600)

hist(plot.data,breaks=100, main=i,xlab="MFI",freq = F,
     ylim=c(0,max(gauss1a,gauss2a)))

if(is.null(neg_label)==F) abline(v=neg.mean,col="springgreen4",lwd=2)
#abline(v=neg.data,col=alpha("springgreen4",0.5),lty=3,lwd=0.2)
if(is.null(blank_label)==F) abline(v=blank.mean,col="darkgrey",lwd=1)
#abline(v=blank.data,col=alpha("darkgrey",0.5),lty=3,lwd=0.2)
abline(v=cut,col="red",lwd=2)
lines(plot_x,gauss1a,col="dodgerblue4")
lines(plot_x,gauss2a,col="dodgerblue4")

legend("topleft",paste0("cutoff (MFI): ",round(cut,3)),lty=1,col="red",cex=0.9,bty="n",y.intersp=0.2,x.intersp=0.2,seg.len=0.5,text.col="red")
  
graphics.off()


