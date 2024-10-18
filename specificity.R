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
library(introdataviz)
library("NatParksPalettes")
## Define a "not" operators -----
`%!in%` <- Negate(`%in%`)
`%!like%` <- Negate(`%like%`)

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
wusmPheno <- read.csv("wusm.pheno.csv", row.names = 1)
humtPheno <- read.csv("humt.pheno.csv",row.names = 1)
humtDLBpheno <- humtPheno[humtPheno$DiseaseStatus!="PD",]
humtPheno <- humtPheno[humtPheno$DiseaseStatus!="DLB",]
humtPheno$DiseaseStatus[humtPheno$DiseaseStatus=="CO"]<-0
humtPheno$DiseaseStatus[humtPheno$DiseaseStatus!=0]<-1
humtPheno$SampleID<-rownames(humtPheno)
humtDLBpheno$DiseaseStatus[humtDLBpheno$DiseaseStatus=="CO"]<-0
humtDLBpheno$DiseaseStatus[humtDLBpheno$DiseaseStatus!=0]<-1
humtDLBpheno$SampleID<-rownames(humtDLBpheno)

#Ensure that each patient has both
#counts and phenotype available
humtPheno<-humtPheno[humtPheno$SampleID%in%rownames(humtCount),]
humtDLBpheno<-humtDLBpheno[humtDLBpheno$SampleID%in%rownames(humtDLBcount),]
humtCount<-humtCount[rownames(humtCount)%in%rownames(humtPheno),]
humtDLBcount<-humtDLBcount[rownames(humtDLBcount)%in%rownames(humtDLBpheno),]
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

# Subset count matrices to transcripts of interest
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

# Set up the training response variable
trainY<-as.numeric(humtPheno$status[sort(rownames(humtPheno)%in%rownames(humtCount))])

# 4 Healthy controls vs AD/DLB/FTD ----
# Test weather the model can differentiate between healthy controls and AD/DLB/FTD
aucdf<-data.frame()
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
  
  trainCount = humtCount[sort(rownames(humtCount)),colnames(humtCount)%in%genestodoml]
  ### AD -----
  # Load data
  load("AD.rlog.sex.age.adj.noFreezerGenes.RData")
  cnt = as.data.frame(t(assay(rld)))
  pheno <- read.csv("AD.DLB.FTD.pheno.csv")
  pheno$pooID <- as.numeric(pheno$Pool.ID) + 106695
  rownames(pheno)<-paste(pheno$pooID,pheno$Library.Barcode,sep = "-")
  pheno<-pheno[rownames(pheno)%in%rownames(cnt),]
  
  ### AD cdr 0.5 k99
  # Subset AD CDR0.5 and controls
  pheno05<-pheno[pheno$Group%in%c("CDR05","Control"),]
  pheno05$Status<-0
  pheno05$Status[pheno05$Group=="CDR05"]<-1
  test<-cnt[rownames(pheno05),]
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  y_test <- as.numeric(pheno05$Status)
  # Testing predictions
  cdr05k99preds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  cdr05k99 = data.frame(Status=y_test,Predictor = as.numeric(cdr05k99preds))
  controlsforspe<-cdr05k99[cdr05k99$Status==0,]
  
  ### AD cdr1 k99
  # Subset AD CDR1 and controls
  pheno1<-pheno[pheno$Group%in%c("CDR1","Control"),]
  pheno1$Status<-0
  pheno1$Status[pheno05$Group=="CDR1"]<-1
  test<-cnt[rownames(pheno1),]
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  y_test <- as.numeric(pheno1$Status)
  # Testing predictions
  cdr1k99preds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  cdr1k99 = data.frame(Status=y_test,Predictor = as.numeric(cdr1k99preds))
  cdr1k99 = cdr1k99[cdr1k99$Status == 1,]
  # Combine all AD
  df<-rbind(cdr05k99,cdr1k99)
  confu = caret::confusionMatrix(data = factor(as.numeric(df$Predictor>0.5)), reference = factor(df$Status))
  auc<-ci.auc(as.numeric(df$Status), as.numeric(df$Predictor),conf.level = 0.9)
  t<-data.frame(Reference="Healthy control",Model=paste(length(genestodoml),"transcripts",sep=" "),
                Status="Alzheimer's disease",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),  
                AUC = paste(round(auc[2],3)," (",round(auc[1],3),"-",round(auc[3],3),")",sep=""),
                PPV=round(confu$byClass[3][[1]],3),NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf,t)
  
  ### DLB -----
  # Load the data
  test = readRDS("DLB.normalized.counts")
  test = as.data.frame(test)
  test = test[test$Status == "DLB",]
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Scale training and testing datasets
  rows<-rownames(x_test)
  x_test <- scaletwosets(trainCount,x_test)
  row.names(x_test)<-rows
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  numeros = as.numeric(test$Status)
  y_test <-replace(numeros, numeros == 2,1)
  # Testing predictions
  dlbpreds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  dlb = data.frame(Status=y_test,Predictor = as.numeric(dlbpreds))
  dlb = dlb[dlb$Status == 1,]
  dlb = rbind(dlb,controlsforspe)
  # Select features
  x_test<-humtDLBcount[sort(rownames(humtDLBcount)),colnames(humtDLBcount)%in%genestodoml]
  # Scale training and testing datasets
  rows<-rownames(x_test)
  x_test <- scaletwosets(trainCount,x_test)
  row.names(x_test)<-rows
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  y_test <-as.numeric(humtDLBpheno$DiseaseStatus)
  # Testing predictions
  dlbpreds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  df = data.frame(Status=y_test, Predictor = as.numeric(dlbpreds))
  df<-rbind(dlb,df)
  confu = caret::confusionMatrix(data = factor(as.numeric(df$Predictor>0.5)), reference = factor(df$Status))
  auc<-ci.auc(as.numeric(df$Status), as.numeric(df$Predictor),conf.level = 0.9)
  t<-data.frame(Reference="Healthy control",Model=paste(length(genestodoml),"transcripts",sep=" "),
                Status="Dementia with Lewy bodies",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),  
                AUC = paste(round(auc[2],3)," (",round(auc[1],3),"-",round(auc[3],3),")",sep=""),
                PPV=round(confu$byClass[3][[1]],3),NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf,t)
  
  ### FTD -----
  test = readRDS("FTD.normalized.counts")
  test = as.data.frame(test)
  test = subset(test, test$Status == "FTD")
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Scale training and testing datasets
  rows<-rownames(x_test)
  x_test <- scaletwosets(trainCount,x_test)
  row.names(x_test)<-rows
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  numeros = as.numeric(test$Status)
  y_test <-replace(numeros, numeros == 2,1)
  # Testing predictions
  ftdpreds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  ftd = data.frame(Status=y_test,Predictor = as.numeric(ftdpreds))
  ftd = ftd[ftd$Status == 1,]
  df = rbind(ftd,controlsforspe)
  confu = caret::confusionMatrix(data = factor(as.numeric(df$Predictor>0.5)), reference = factor(df$Status))
  auc<-ci.auc(as.numeric(df$Status), as.numeric(df$Predictor),conf.level = 0.9)
  t<-data.frame(Reference="Healthy control",Model=paste(length(genestodoml),"transcripts",sep=" "),
                Status="Frontotemporal dementia",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),  
                AUC = paste(round(auc[2],3)," (",round(auc[1],3),"-",round(auc[3],3),")",sep=""),
                PPV=round(confu$byClass[3][[1]],3),NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf,t)
  
}

# 5 PD vs AD/DLB/FTD -----
# Test weather the model can differentiate between PD and AD/DLB/FTD
aucdf<-data.frame()
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
  
  # Set up testing data and subset only genes in dfkld
  testX <- wusmCount
  testX = testX[sort(rownames(testX)),colnames(testX)%in%genestodoml]
  testX <- as.matrix(testX)
  testY<-as.numeric(wusmPheno$status[sort(rownames(wusmPheno)%in%rownames(wusmCount))])
  # Testing preditions
  predsTest<-predict(bestRidge, newx = testX, type = "response", standardize = F, intercept = F)
  df<-data.frame(Status=testY, Predictor=as.numeric(predsTest))
  PDpreds<-df[df$Status==1,]
  trainCount = humtCount[sort(rownames(humtCount)),colnames(humtCount)%in%genestodoml]
  ### AD ----
  # Load data
  load("AD.rlog.sex.age.adj.noFreezerGenes.RData")
  cnt = as.data.frame(t(assay(rld)))
  pheno <- read.csv("AD.DLB.FTD.pheno.csv")
  pheno$pooID <- as.numeric(pheno$Pool.ID) + 106695
  rownames(pheno)<-paste(pheno$pooID,pheno$Library.Barcode,sep = "-")
  pheno<-pheno[rownames(pheno)%in%rownames(cnt),]
  
  ### AD cdr 0.5 k99
  # Subset AD CDR0.5 and controls
  pheno05<-pheno[pheno$Group%in%c("CDR05","Control"),]
  pheno05$Status<-0
  pheno05$Status[pheno05$Group=="CDR05"]<-1
  test<-cnt[rownames(pheno05),]
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  y_test <- as.numeric(pheno05$Status)
  # Testing predictions
  cdr05k99preds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  cdr05k99 = data.frame(Status=y_test,Predictor = as.numeric(cdr05k99preds))
  
  ### AD cdr1 k99
  # Subset AD CDR1 and controls
  pheno1<-pheno[pheno$Group%in%c("CDR1","Control"),]
  pheno1$Status<-0
  pheno1$Status[pheno1$Group=="CDR1"]<-1
  test<-cnt[rownames(pheno1),]
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  y_test <- as.numeric(pheno1$Status)
  # Testing predictions
  cdr1k99preds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  cdr1k99 = data.frame(Status=y_test,Predictor = as.numeric(cdr1k99preds))
  # Combine all AD
  df<-rbind(cdr05k99,cdr1k99)
  df<-df[df$Status==1,]
  df$Status<-0
  df<-rbind(df,PDpreds)
  confu = caret::confusionMatrix(data = factor(as.numeric(df$Predictor>0.5)), reference = factor(df$Status))
  auc<-ci.auc(as.numeric(df$Status), as.numeric(df$Predictor),conf.level = 0.9)
  t<-data.frame(Reference="Parkinson's disease",Model=paste(length(genestodoml),"transcripts",sep=" "),
                Status="Alzheimer's disease",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),  
                AUC = paste(round(auc[2],3)," (",round(auc[1],3),"-",round(auc[3],3),")",sep=""),
                PPV=round(confu$byClass[3][[1]],3),NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf,t)
 
  ### DLB -----
  # Load the data
  test = readRDS("DLB.normalized.counts")
  test = as.data.frame(test)
  test = test[test$Status == "DLB",]
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Scale training and testing datasets
  rows<-rownames(x_test)
  x_test <- scaletwosets(trainCount,x_test)
  row.names(x_test)<-rows
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  numeros = as.numeric(test$Status)
  y_test <-replace(numeros, numeros == 2,1)
  # Testing predictions
  dlbpreds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  dlb = data.frame(Status=y_test,Predictor = as.numeric(dlbpreds))
  dlb = dlb[dlb$Status == 1,]
  # Select features
  x_test<-humtDLBcount[sort(rownames(humtDLBcount)),colnames(humtDLBcount)%in%genestodoml]
  # Scale training and testing datasets
  rows<-rownames(x_test)
  x_test <- scaletwosets(trainCount,x_test)
  row.names(x_test)<-rows
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  y_test <-as.numeric(humtDLBpheno$DiseaseStatus)
  # Testing predictions
  dlbpreds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  df = data.frame(Status=y_test, Predictor = as.numeric(dlbpreds))
  df<-rbind(dlb,df)
  df<-df[df$Status==1,]
  df$Status<-0
  df<-rbind(df,PDpreds)
  confu = caret::confusionMatrix(data = factor(as.numeric(df$Predictor>0.5)), reference = factor(df$Status))
  auc<-ci.auc(as.numeric(df$Status), as.numeric(df$Predictor),conf.level = 0.9)
  t<-data.frame(Reference="Parkinson's disease",Model=paste(length(genestodoml),"transcripts",sep=" "),
                Status="Dementia with Lewy bodies",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),  
                AUC = paste(round(auc[2],3)," (",round(auc[1],3),"-",round(auc[3],3),")",sep=""),
                PPV=round(confu$byClass[3][[1]],3),NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf,t)
  
  ### FTD -----
  test = readRDS("FTD.normalized.counts")
  test = as.data.frame(test)
  test = subset(test, test$Status == "FTD")
  # Select features
  x_test<-test[,colnames(test)%in%genestodoml]
  # Scale training and testing datasets
  rows<-rownames(x_test)
  x_test <- scaletwosets(trainCount,x_test)
  row.names(x_test)<-rows
  # Compute Z score
  x_test <- as.data.frame(lapply(x_test, function(x) (x - mean(x))/sd(x) ))
  x_test <- as.matrix(x_test)
  # Code status
  numeros = as.numeric(test$Status)
  y_test <-replace(numeros, numeros == 2,1)
  # Testing predictions
  ftdpreds <- predict(bestRidge, newx = x_test, type = "response" , standarize = F, intercept = F)
  # Model roc comes from above code ( genes model)
  ftd = data.frame(Status=y_test,Predictor = as.numeric(ftdpreds))
  df = ftd[ftd$Status == 1,]
  df$Status<-0
  df<-rbind(df,PDpreds)
  confu = caret::confusionMatrix(data = factor(as.numeric(df$Predictor>0.5)), reference = factor(df$Status))
  auc<-ci.auc(as.numeric(df$Status), as.numeric(df$Predictor),conf.level = 0.9)
  t<-data.frame(Reference="Parkinson's disease",Model=paste(length(genestodoml),"transcripts",sep=" "),
                Status="Frontotemporal dementia",BestAccuracy=round(confu$byClass[11][[1]],3),
                CohenKappa=round(confu$overall[2][[1]],3),Sensitivity=round(confu$byClass[1][[1]],3),
                Specificity=round(confu$byClass[2][[1]],3),  
                AUC = paste(round(auc[2],3)," (",round(auc[1],3),"-",round(auc[3],3),")",sep=""),
                PPV=round(confu$byClass[3][[1]],3),NPV=round(confu$byClass[4][[1]],3))
  aucdf<-rbind(aucdf,t)
  
}
