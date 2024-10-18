#!/usr/bin/env Rscript
# 0 Prepare the environment -----
## Clear the environment -----
rm(list = ls())
lapply(paste0("package:", names(sessionInfo()$otherPkgs)), detach, character.only = TRUE, unload = TRUE)
dev.off()
## Load the necessary libraries -----
library(arules)
library(entropy)
library(dplyr)
library(scales)
library(data.table)
library(glmnet)
library(ggplot2)
library(pROC)
library("NatParksPalettes")
# Define a "not in" opperator -----
`%!in%` <- Negate(`%in%`)

# 1 Define a function to scale two different datasets -----
k = 1
scaletwosets = function(reference, dftoscale){
  # As data frame
  reference = as.data.frame(reference)
  dftoscale = as.data.frame(dftoscale)
  reference = reference[,colnames(reference) %in% colnames(dftoscale)]
  # Make data frame with the same number of rows and columns
  dfscaled = data.frame(matrix(0, nrow = nrow(dftoscale), ncol = ncol(dftoscale),
                                 dimnames = list(NULL, colnames(dftoscale))) )
  for(k in 1:length(colnames(reference))){
    # Compute max and min for each feature in the training reference
    df1range = NULL
    df1range = range(reference[,k])
    df1rangemin = df1range[1]
    df1rangemax = df1range[2]
    # Scale to this min and max the values in the dftoscale
    vec_range <- rescale(dftoscale[,k], to = c(df1rangemin, df1rangemax))
    dfscaled[,k] = vec_range
  }
  return(dfscaled)
}
# 2 Define a function to compute the Kullback–Leibler divergence (KL) -----
KLgenedf <- function(df1, df2, numBins = 5)
{
  # Order an obtain the same genes  
  genes <- colnames(df1)
  df1 <- df1[, genes]
  df2 <- df2[, genes]
  # Find the range of values (min and max) for each dataset
  dfkl = data.frame(matrix(0, nrow = 1, ncol = length(genes), dimnames = list(NULL, genes)) )
  for(gene in genes) {
    df1range = range(df1[,gene])
    df2range = range(df2[,gene])
    minrange = min(df1range[1], df2range[1])
    maxrange = max(df1range[2], df2range[2])
    ourrange = c(minrange,maxrange)
    # Convert continouus variable that is gene distribution
    # into a categorical variable (factor)
    disdf1 = discretize(df1[,gene], numBins = numBins, r = ourrange)
    disdf2 = discretize(df2[,gene], numBins = numBins, r = ourrange)
    
    disdf1 = disdf1/sum(disdf1)
    disdf2 = disdf2/sum(disdf2)
    # Plug-In Estimator of the Kullback-Leibler 
    # divergence and of the Chi-Squared Divergence
    kldf1 = KL.plugin(disdf1,disdf2)
    kldf2 = KL.plugin(disdf2, disdf1)
    kldif = kldf1+kldf2
    
    dfkl[,gene] <- kldif
    
  }
  return(dfkl)
}

# 3 Load and format the data -----
## Load meta analysis results -----
df<-read.csv("meta.humt.wusm.DE.csv")
genes<-df$Gene
df<-as.data.frame(sapply(df[,-1],as.numeric))
df$Gene<-genes
DEgenes<-df$Gene[!is.na(df$MetaPadj)&df$MetaPadj<0.05]

## Load rlog transformed and normalized count data -----
load("wusm.rlog.sex.age.adj.noFreezerGenes.noMedGenes.RData")
wusmCount<-as.data.frame(t(assay(rld)))
load("humt.rlog.sex.age.adj.noFreezerGenes.noMedGenes.RData")
humtCount<-as.data.frame(t(assay(rld)))
## Load phenotype data -----
wusmPheno <- read.csv("wusm.pheno.csv", row.names = 1)
humtPheno <- read.csv("humt.pheno.csv",row.names = 1)
humtPheno <- humtPheno[humtPheno$DiseaseStatus!="DLB",]
humtPheno$DiseaseStatus[humtPheno$DiseaseStatus=="CO"]<-0
humtPheno$DiseaseStatus[humtPheno$DiseaseStatus!=0]<-1
humtPheno$SampleID<-rownames(humtPheno)
# Ensure that each patient has both
# counts and phenotype available
humtPheno<-humtPheno[humtPheno$SampleID%in%rownames(humtCount),]
humtCount<-humtCount[rownames(humtCount)%in%rownames(humtPheno),]
wusmPheno<-wusmPheno[rownames(wusmPheno)%in%rownames(wusmCount),]
wusmCount<-wusmCount[rownames(wusmCount)%in%rownames(wusmPheno),]
# Subset each dataset to contain only genes that 
# have the same sign in beta between the two datasets
wusmCount<-wusmCount[,colnames(wusmCount)%in%DEgenes]
humtCount<-humtCount[,colnames(humtCount)%in%DEgenes]
# Scale discvery and replication datasets
rows<-rownames(wusmCount)
wusmCount <- scaletwosets(humtCount,wusmCount)
row.names(wusmCount)<-rows
# 4 Calculate Z scores and Kullback-Leibler divergence -----
rows<-rownames(wusmCount)
wusmCount <- as.data.frame(lapply(wusmCount, function(x) (x - mean(x))/sd(x) ))
row.names(wusmCount)<-rows
rows<-rownames(humtCount)
humtCount <- as.data.frame(lapply(humtCount, function(x) (x - mean(x))/sd(x) ))
row.names(humtCount)<-rows
kldiv = KLgenedf(wusmCount, humtCount)
tkldiv<-as.data.frame(t(kldiv))
names(tkldiv)<-"kld"
# 5 Select genes based on Kullback-Leibler divergence -----
klsdf= seq(0.01, 1, by= 0.01)
bestacu<-data.frame()

for(i in  1:length(klsdf))
{
  klsd =  klsdf[i]
  dfkld = kldiv[,kldiv < klsd]
  genestodoml = colnames(dfkld)
  #bestacu[i,3] = length(genestodoml)
  if(length(genestodoml!=0))
  {
    # Set up training data and subset only genes in dfkld
    trainX <- humtCount
    trainX = trainX[sort(rownames(trainX)),colnames(trainX)%in%genestodoml]
    trainX <- as.matrix(trainX)
    # Code status
    trainY<-as.numeric(humtPheno$DiseaseStatus[sort(rownames(humtPheno)%in%rownames(humtCount))])
    # Set seed
    set.seed(567)
    # Train ridge
    cvRidge <- cv.glmnet(x = trainX, y = factor(trainY), alpha = 0,
                          nfolds = 5, type.measure = "mse", standardize = F,
                          family = "binomial", intercept = F)
    # Best model lambda
    bestRidge <- glmnet(x = trainX, y = trainY, alpha = 0, lambda = cvRidge$lambda.min,
                     standardize = F, family = "binomial", intercept = F)
    # Training predictions
    predsTrain<-predict(bestRidge, newx = trainX, type = "response", standardize = F)
    # MSE (Training)
    mseTrain <- mean((predsTrain - trainY)^2)
    
    confu = caret::confusionMatrix(data = factor(as.numeric(predsTrain>0.5)), reference = factor(trainY))
    auc<-ci.auc(as.numeric(trainY), as.numeric(predsTrain),conf.level = 0.9)
    t<-data.frame(KLDvalue=klsd,nGenes=length(genestodoml),Status="Training",BestAccuracy=round(confu$byClass[11][[1]],3),
                  CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                  Specificity=round(confu$byClass[2][[1]],3),AUClow = round(auc[1],3), 
                  AUC = round(auc[2],3), AUChigh = round(auc[3],3),PPV=round(confu$byClass[3][[1]],3),
                  NPV=round(confu$byClass[4][[1]],3))
    
    bestacu = rbind(bestacu,t)
  }

}

# Selected KLD cutoff values: 0.1, 0.23 and 0.61