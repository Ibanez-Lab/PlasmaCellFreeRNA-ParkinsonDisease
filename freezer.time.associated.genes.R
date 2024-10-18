#!/usr/bin/env Rscript
# 0 Prepare the environment -----
## Load the necessary libraries -----
library(dplyr)
library(DESeq2)
library(reshape2)
library(data.table)

## Define a "not in" opperator -----
`%!in%` <- Negate(`%in%`)

# 1 Load and format the data -----
# Load the count matrix
countMatrix <- read.csv("cfRNA.count.matrix.csv",row.names = 1,check.names = F)
# Load the phenotype information
pheno <- read.csv("pheno.csv")

# Remove any genes with less than 10 counts in 90% or more of samples
countMatrix <- countMatrix[(matrixStats::rowCounts(countMatrix<10) < round(0.9*dim(countMatrix)[2])),]
# Remove outliers
countMatrix <- countMatrix[,colnames(countMatrix) %!in% outliers]

# Subset only control samples
pheno<-pheno[pheno$DiseaseStatus=="CO",]

#Ensure that phenotype is available for all samples
pheno<-pheno[pheno$SampleID%in%colnames(countMatrix),]
countMatrix<-countMatrix[,colnames(countMatrix)%in%pheno$SampleID]
pheno$Pool<-as.character(pheno$Pool)
row.names(pheno)<-pheno$SampleID

# Ensure that count and phenotype data are in the same order
pheno<-pheno[colnames(countMatrix),]

#Ensure count matrix is integer
rows<-row.names(countMatrix)
countMatrix<-sapply(countMatrix, as.integer)
row.names(countMatrix)<-rows

# 2 Find differentially expressed genes -----
DDS<-DESeqDataSetFromMatrix(countData = countMatrix, colData = pheno, design = ~ AgeDraw+factor(Sex)+FreezerTime_years)
DE <-results(DESeq(DDS), tidy = T, independentFiltering=FALSE)
# Save DE results
write.csv(DE, "DE.age.sex.freezerTime.csv", row.names = F)
