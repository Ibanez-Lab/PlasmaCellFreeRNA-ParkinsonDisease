#!/usr/bin/env Rscript
# 0 Prepare the environment -----
## Load the necessary libraries -----
library(dplyr)
library(DESeq2)
library(reshape2)
library(data.table)
library(flexiblas)
n <- flexiblas_load_backend("OPENBLAS-SERIAL")
flexiblas_switch(n)
## Define a "not in" opperator -----
`%!in%` <- Negate(`%in%`)

# 1 Load and format the data -----
## Load the necessary data -----
# Load the count matrix
countMatrix <- read.csv("cfRNA.count.matrix.csv",row.names = 1,check.names = F)
# Load the phenotype information
pheno <- read.csv("pheno.csv")
# Load freezer time associated genes
freezergenes<-read.csv("DE.age.sex.freezerTime.csv")
freezergenes<-freezergenes$row[!is.na(freezergenes$pvalue)&freezergenes$pvalue<0.05]
# Load PD medication associated genes
medgenes<-read.csv("DE.CAonly.age.sex.PDmeds.csv")
medgenes<-medgenes$row[!is.na(medgenes$pvalue)&medgenes$pvalue<0.05]

## Format the data -----
# Remove genes associated with freezer time or medication usage
countMatrix<-countMatrix[rownames(countMatrix)%!in%c(freezergenes,medgenes),]
# Remove any genes with less than 10 counts in 90% or more of samples
countMatrix <- countMatrix[(matrixStats::rowCounts(countMatrix<10) < round(0.9*dim(countMatrix)[2])),] #removed 32301 genes
# Remove outliers
countMatrix <- countMatrix[,colnames(countMatrix) %!in% outliers]

# Change DiseaseStatus to a binary (0/1) variable
pheno$DiseaseStatus[pheno$DiseaseStatus=="CO"] <- 0
pheno$DiseaseStatus[pheno$DiseaseStatus!=0] <- 1

#Ensure that phenotype is available for all samples
pheno<-pheno[pheno$SampleID%in%colnames(countMatrix),]
rownames(pheno)<-pheno$SampleID
countMatrix<-countMatrix[,colnames(countMatrix)%in%pheno$SampleID]
pheno<-pheno[colnames(countMatrix),]

# Ensure that values in the count matrix are integers
rows<-row.names(countMatrix)
countMatrix<-sapply(countMatrix, as.integer)
row.names(countMatrix)<-rows
pheno$Pool<-as.character(pheno$Pool)

# 2 Differential Expression Analysis -----
DDS<-DESeqDataSetFromMatrix(countData = countMatrix, colData = pheno, design = ~ AgeDraw+factor(Sex)+factor(DiseaseStatus))
DE <-results(DESeq(DDS), tidy = T, independentFiltering=FALSE)
# Save DE results
write.csv(DE, "DE.noMedGenes.noFreezerGenes.age.sex.status.csv", row.names = F)
