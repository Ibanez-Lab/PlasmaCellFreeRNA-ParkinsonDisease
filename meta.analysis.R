#!/usr/bin/env Rscript
# 0 Prepare the environment -----
## Load the necessary libraries -----
library(metaRNASeq)
library(ggplot2)
library("NatParksPalettes")
# 1 Load and format the data -----
humt<-read.csv("humt.DE.noMedGenes.noFreezerGenes.age.sex.status.csv")
wusm<-read.csv("wusm.DE.noMedGenes.noFreezerGenes.age.sex.status.csv")
# Format and subset the data
humt<-humt[,c(1,3,6)]
wusm<-wusm[,c(1,3,6)]
names(humt)<-c("Gene", "humtFC", "humtPraw")
names(wusm)<-c("Gene", "wusmFC", "wusmPraw")
DE<-merge(humt, wusm, by="Gene")
DE<-DE[!is.na(DE$humtPraw),]
DE<-DE[!is.na(DE$wusmPraw),]
nonDE<-DE[sign(DE$humtFC)!=sign(DE$wusmFC),]
DE<-DE[sign(DE$humtFC)==sign(DE$wusmFC),]
# 2 Perform meta-analyses -----
fishcomb <- fishercomb(list(DE$humtPraw,DE$wusmPraw), BHth = 0.05)
# Format results
metaDE<-cbind(DE,fishcomb$rawpval,fishcomb$adjpval)
names(metaDE)[c(6,7)]<-c("MetaPraw", "MetaPadj")
nrow(metaDE[metaDE$MetaPadj<0.05,])
nonDE$MetaPraw<-"-"
nonDE$MetaPadj<-"-"
df<-rbind(metaDE,nonDE)
write.csv(df, "meta.humt.wusm.DE.csv", row.names = F)
