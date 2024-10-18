#!/usr/bin/env Rscript
# 0 Prepare the environment -----
## Load the necessary libraries ------
# library(tximport)
library(DESeq2)
library(dplyr)
library(data.table)
library(flexiblas)
n <- flexiblas_load_backend("OPENBLAS-SERIAL")
flexiblas_switch(n)
## Define a "not in" opperator -----
`%!in%` <- Negate(`%in%`)

# 1 Load and format the data -----
# Load the count matrix 
countMatrix <- read.csv("cfRNA.count.matrix.csv",row.names = 1,check.names = F)
#Load freezer time associated genes
freezergenes<-read.csv("DE.age.sex.FreezerTime.csv")
freezergenes<-freezergenes$row[freezergenes$row<0.05]
#Load phenotype information for the samples
pheno <- read.csv("pheno.csv", check.names = F, row.names = 1)
#Summarize PD medication info
pheno<-pheno[pheno$DiseaseStatus=="PD",]
pheno$PDmed<-0
pheno$PDmed[pheno$LevoDopa=="Y"|pheno$DopaminePromoter=="Y"]<-1

# Remove outliers
countMatrix<-countMatrix[,colnames(countMatrix)%!in%outliers,]
# Remove freezder time associated genes
countMatrix<-countMatrix[rownames(countMatrix)%!in%freezergenes,]
#Ensure that phenotype data is available for all samples
countMatrix <- countMatrix[,colnames(countMatrix)%in%rownames(pheno)]
pheno<-pheno[rownames(pheno)%in%colnames(countMatrix),]
pheno<-pheno[colnames(countMatrix),]

# 2 Find differentially expressed genes
DDS <- DESeqDataSetFromMatrix(countData = countMatrix, colData = pheno, design = ~ age_at_draw+factor(sex)+factor(PDmed))
DE <-results(DESeq(DDS), tidy = T, independentFiltering=FALSE)
write.csv(DE,"DE.CAonly.age.sex.PDmeds.csv", row.names=F)
