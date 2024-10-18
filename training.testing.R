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
library(pROC)
library(ggplot2)
library("NatParksPalettes")
## Define a "not in" operator -----
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

# 2 Load and format the data -----
## Load rlog transformed and normalized count data -----
load("wusm.rlog.sex.age.adj.noFreezerGenes.noMedGenes.RData")
wusmCount<-as.data.frame(t(assay(rld)))
load("humt.rlog.sex.age.adj.noFreezerGenes.noMedGenes.RData")
humtCount<-as.data.frame(t(assay(rld)))
humtDLBcount<-humtCount
## Load phenotype data -----
wusmPheno <- read.csv("/03-DryLab/04-Analyses/2020_cfRNA-Biomarkers_LI/2020_PD_cfRNA_ABeric/pheno.medication.SV1-10.smplSVA.ComBatSeq.rm4outliers.csv", row.names = 1)
humtPheno <- read.csv("/03-DryLab/04-Analyses/2020_cfRNA-Biomarkers_LI/2023_Tatooine_ABeric/HUMT_pheno_SVA.csv",row.names = 1)
humtPheno <- humtPheno[humtPheno$DiseaseStatus!="DLB",]
humtPheno$DiseaseStatus[humtPheno$DiseaseStatus=="CO"]<-0
humtPheno$DiseaseStatus[humtPheno$DiseaseStatus!=0]<-1
humtPheno$SampleID<-rownames(humtPheno)
#Ensure that each patient has both
#counts and phenotype available
humtPheno<-humtPheno[humtPheno$SampleID%in%rownames(humtCount),]
humtCount<-humtCount[rownames(humtCount)%in%rownames(humtPheno),]
wusmPheno<-wusmPheno[rownames(wusmPheno)%in%rownames(wusmCount),]
wusmCount<-wusmCount[rownames(wusmCount)%in%rownames(wusmPheno),] 
# Subset and format phenotype data
humtPheno<-humtPheno[,c(2,6,7)]
wusmPheno<-wusmPheno[,c(11,12,10)]
names(humtPheno)<-c("age","sex","status")
names(wusmPheno)<-c("age","sex","status")
## Load feature selection results -----
geneKLD<-read.csv("KLD.value.DE.transcripts.csv")
names(geneKLD)[1]<-"Gene"
lims<-c(0.1,0.23,0.61) # number of genes: 26, 87 and 191, respectively

# Subset count matrices to genes of interest
humtCount<-humtCount[,colnames(humtCount)%in%geneKLD$Gene]
wusmCount<-wusmCount[,colnames(wusmCount)%in%geneKLD$Gene]

# 3 Scale discovery and replication datasets -----
rows<-rownames(wusmCount)
wusmCount <- scaletwosets(humtCount,wusmCount)
row.names(wusmCount)<-rows
# Calculate Z scores
rows<-rownames(wusmCount)
wusmCount <- as.data.frame(lapply(wusmCount, function(x) (x - mean(x))/sd(x) ))
row.names(wusmCount)<-rows
rows<-rownames(humtCount)
humtCount <- as.data.frame(lapply(humtCount, function(x) (x - mean(x))/sd(x) ))
row.names(humtCount)<-rows

# 4 Calculate AUC for age and sex as baseline -----
allPDpheno<-rbind(wusmPheno,humtPheno)
fit <- glm(factor(status) ~ age+factor(sex), data = allPDpheno,family=binomial())
preds = predict(fit, type="response")
auc<-ci.auc(as.numeric(allPDpheno$status), as.numeric(preds),conf.level = 0.9)
aucdf<-data.frame(Genes="age+sex",Status="Baseline",BestAccuracy=NA,CohenKappa=NA,Sensitivity=NA,
                  Specificity=NA,AUClow = round(auc[1],3),AUC = round(auc[2],3), AUChigh = round(auc[3],3),
                  PPV=NA,NPV=NA)

# 5 Test predictive models with the selected transcripts -----

# Set up training an testing response variables
trainY<-as.numeric(humtPheno$status[sort(rownames(humtPheno)%in%rownames(humtCount))])
testY<-as.numeric(wusmPheno$status[sort(rownames(wusmPheno)%in%rownames(wusmCount))])

for (lim in lims)
{
  # Select genes within the desired kld cutoff
  genestodoml<-geneKLD$Gene[geneKLD$kld<lim]
  # Set up training data and subset only genes of interest
  trainX <- humtCount
  trainX = trainX[sort(rownames(trainX)),colnames(trainX)%in%genestodoml]
  trainX <- as.matrix(trainX)
  # Set seed
  set.seed(567)
  # Train ridge
  cvRidge <- cv.glmnet(x = trainX, y = factor(trainY), alpha = 0,
                       nfolds = 5, type.measure = "mse", standardize = F,
                       family = "binomial", intercept = F)
  # Best model lambda
  bestRidge <- glmnet(x = trainX, y = trainY, alpha = 0, lambda = cvRidge$lambda.min,
                      standardize = F, family = "binomial", intercept = F)
  # Compute coeficients
  dfCoef <- coef(bestRidge) %>% as.matrix() %>% as_tibble(rownames = "predictor")
  # Training predictions
  predsTrain<-predict(bestRidge, newx = trainX, type = "response", standardize = F)
  # Confusion Matrix
  confu = caret::confusionMatrix(data = factor(as.numeric(predsTrain>0.5)), reference = factor(trainY))
  auc<-ci.auc(as.numeric(trainY), as.numeric(predsTrain),conf.level = 0.9)
  t<-data.frame(Genes=length(genestodoml),Status="Training",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),AUClow = round(auc[1],3), 
                AUC = round(auc[2],3), AUChigh = round(auc[3],3),PPV=round(confu$byClass[3][[1]],3),
                NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf,t)
  # Set up testing data and subset only genes in dfkld
  testX <- wusmCount
  testX = testX[sort(rownames(testX)),colnames(testX)%in%genestodoml]
  testX <- as.matrix(testX)
  # Testing preditions
  predsTest<-predict(bestRidge, newx = testX, type = "response", standardize = F, intercept = F)
  # Confusion Matrix
  confu = caret::confusionMatrix(data = factor(as.numeric(predsTest>0.5)), reference = factor(testY))
  auc<-ci.auc(as.numeric(testY), as.numeric(predsTest),conf.level = 0.9)
  t<-data.frame(Genes=length(genestodoml),Status="Testing",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),AUClow = round(auc[1],3), 
                AUC = round(auc[2],3), AUChigh = round(auc[3],3),PPV=round(confu$byClass[3][[1]],3),
                NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf, t)
}
