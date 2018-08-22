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

fit.ab2<-normalmixEM(antibody1, lambda=c(0.5,0.5),mu=c(4.0,6.0),k=2)

#5 Plot Cut off Value
summary(fit.ab2)

cutoff<-uniroot(f,c(0,15),prob=0.90,lambda=fit.ab2$lambda,mu=fit.ab2$mu,sigma=fit.ab2$sigma,k=2,k1=1)$root

cutoff2<-uniroot(f2,c(0,15),prob=0.90,lambda=fit.ab2$lambda,mu=fit.ab2$mu,sigma=fit.ab2$sigma,k=2,k1=2)$root

sum(antibody1>cutoff2)/length(antibody1)

sum(antibody1<cutoff)/length(antibody1)

1-sum(antibody1>cutoff2)/length(antibody1)-sum(antibody1<cutoff)/length(antibody1)

plot(antibody1,fit.ab2$posterior[,1],type='n',xlim=c(0,10),lwd=2,ylim=c(0,1),col='green',las=1,xlab='wb titres',ylab='classification probability')

rect(cutoff,-0.04,cutoff2,1.04,col='light grey',lwd=1.5)

title('C',adj=0,cex.main=1.5)

lines(antibody1,fit.ab2$posterior[,2],lwd=2,col='red')

abline(h=0.90,lwd=1.5,lty=2)

abline(v=cutoff,lwd=1)

abline(v=cutoff2,lwd=1)

lines(antibody1,fit.ab2$posterior[,1],lwd=2,col='blue')

lines(antibody1,fit.ab2$posterior[,2],lwd=2,col='purple')

legend(1125,0.8,c(expression(S^'-'),expression(S^'+')),lty=c(1,1),lwd=2,col=c('blue','purple'))

