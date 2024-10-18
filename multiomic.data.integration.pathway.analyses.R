#!/usr/bin/env Rscript
# 0 Prepare the environment -----
## Clear the environment -----
rm(list = ls())
lapply(paste0("package:", names(sessionInfo()$otherPkgs)), detach, character.only = TRUE, unload = TRUE)
dev.off()
## Load the necessary libraries -----
library(metaRNASeq)
library(readxl)
library(ggplot2)
library("NatParksPalettes")
library(data.table)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
## Define "not ()" opperators -----
`%!in%` <- Negate(`%in%`) 
`%!like%` <- Negate(`%like%`)
# cfRNA -----
cfrna<-read.csv("meta.humt.wusm.DE.csv")
cfrna$Gene<-sapply(cfrna$Gene, function(a) strsplit(a,"\\.")[[1]][1])
cfrna$MetaPadj<-as.numeric(cfrna$MetaPadj)
hsmart<-readRDS("biomart.v99.ensemblGeneID.entrezGeneID.hgncSymbol.chr.start.end.RData")
hsmart<-getBM(attributes = c('chromosome_name','start_position','end_position','ensembl_gene_id', 'entrezgene_id',
                             'hgnc_symbol'), filters = 'ensembl_gene_id', values = cfrna$Gene, mart = hsmart)
cfrna<-merge(cfrna,hsmart, by.x="Gene", by.y="ensembl_gene_id")

outGO<-enrichGO(cfrna$entrezgene_id[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05], 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
outGO<-outGO[,c(2:4,6)]
outKEGG<-enrichKEGG(cfrna$entrezgene_id[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05], pvalueCutoff=0.05)@result
outKEGG<-outKEGG[,c(3:5,7)]

#Set up variables for hypergeometric test
K<-60721 # number of unique genes in gencode v33
N<-as.numeric(length(unique(cfrna$Gene[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05]))) # number of differentially expressed transcripts

# Find pathways that genes nominally significant in both HUMT and WUSM are enriched in
GOterms<-enrichGO(cfrna$entrezgene_id[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$disPraw<0.05&cfrna$repPraw<0.05], 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(cfrna$entrezgene_id[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$disPraw<0.05&cfrna$repPraw<0.05], pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)

# I Whole blood RNA seq -----
# Load and format the data
blood<-read.csv("whole.blood.DE.metaanalysis.csv")
blood<-blood[sign(blood$disFC)==sign(blood$repFC)&!is.na(blood$Gene),]
# Merge with cfRNA data
df<-merge(cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05,c(12,11,2,7)],
          blood[!is.na(blood$fishPraw)&blood$fishPraw<0.05,c(1,2,8)],by.x="hgnc_symbol",by.y="Gene")
df<-df[sign(df$disFC.x)==sign(df$disFC.y),]
# Perform pathway enrichment analyses
GOterms<-enrichGO(df$entrezgene_id, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(df$entrezgene_id, pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)
# Perform hypergeometric test to evaluate significance of overlap
k<-as.numeric(length(unique(blood$Gene[!is.na(blood$fishPraw)&blood$fishPraw<0.05])))
n<-as.numeric(length(unique(df$hgnc_symbol)))
phyper(n-1, k, K-k, N, lower.tail = F)

# II plasma proteome -----
# load the plasma proteomic data and merge with cfRNA data
plasmaP<-read.csv("plasma.proteomic.PDvsCTRL.csv")
df<-merge(cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05,c(12,11,2,7)],
          plasmaP[!is.na(plasmaP$Pvalue)&plasmaP$Pvalue<0.05,c(9,10,12)],by.x="hgnc_symbol",by.y="Protein")
# Perform hypergeometric test to evaluate significance of overlap
k<-as.numeric(length(plasmaP$Protein[!is.na(plasmaP$Pvalue)&plasmaP$Pvalue<0.05]))
n<-as.numeric(length(df$hgnc_symbol))
phyper(n-1, k, K-k, N, lower.tail = F)

# III brain RNAseq -----
brain<-read.csv("brain.DE.csv", row.names = 1)
df<-merge(cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05,c(1,11,2,7)],
          brain[!is.na(brain$pvalue)&brain$pvalue<0.05,c(1,3,6)],by.x="Gene",by.y="ensembl_gene_id")
# Perform pathway enrichment analyses
GOterms<-enrichGO(df$entrezgene_id, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(df$entrezgene_id, pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)
# Perform hypergeometric test to evaluate significance of overlap
k<-as.numeric(length(unique(brain$ensembl_gene_id[!is.na(brain$pvalue)&brain$pvalue<0.05])))
n<-as.numeric(length(unique(df$Gene)))
phyper(n-1, k, K-k, N, lower.tail = F)
# IV CSF proteome -----
csfP<-read.csv("CSF.proteomic.PDvsCTRL.csv")
df<-merge(cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05,c(12,11,2,7)],
          csfP[!is.na(csfP$Pvalue)&csfP$Pvalue<0.05,c(9,10,12)],by.x="hgnc_symbol",by.y="Protein")
# Perform pathway enrichment analyses
GOterms<-enrichGO(df$entrezgene_id, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(df$entrezgene_id, pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)
# Perform hypergeometric test to evaluate significance of overlap
k<-as.numeric(length(csfP$Protein[!is.na(csfP$Pvalue)&csfP$Pvalue<0.05]))
n<-as.numeric(length(df$hgnc_symbol))
phyper(n-1, k, K-k, N, lower.tail = F)

# V anterior cingulate cortex snRNAseq -----
corticalSC<-read.csv("Felke.2021.anterior.cingulate.cortex.DEG.summary.csv")
df<-cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$hgnc_symbol%in%corticalSC$hgnc_symbol[corticalSC$Disease=="PD"],]
# Perform pathway enrichment analyses
GOterms<-enrichGO(df$entrezgene_id, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(df$entrezgene_id, pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)

# VI subcortical putamen snRNAseq -----
subcorticalSC<-read.csv("Xu.2023.PD.subcortical.putamen.DEG.summary.csv")
df<-cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$hgnc_symbol%in%subcorticalSC$hgnc_symbol[subcorticalSC$Disease=="PD"],]
# Perform pathway enrichment analyses
GOterms<-enrichGO(df$entrezgene_id, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(df$entrezgene_id, pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)

# VII Wang et al. substantia nigra -----
sn1<-read.csv("Wang.2024.PD.substantia.nigra.DEG.summary.csv")
df<-cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$hgnc_symbol%in%sn1$hgnc_symbol,]
# Perform pathway enrichment analyses
GOterms<-enrichGO(df$entrezgene_id, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(df$entrezgene_id, pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)

# VIII in-house substantia nigra -----
sn2<-read.csv("in.house.PD.substantia.nigra.DEG.summary.csv")
df<-cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$hgnc_symbol%in%sn2$hgnc_symbol,]
# Perform pathway enrichment analyses
GOterms<-enrichGO(df$entrezgene_id, 'org.Hs.eg.db', ont="ALL", pvalueCutoff=0.05)@result
KEGGterms<-enrichKEGG(df$entrezgene_id, pvalueCutoff=0.05)@result
outKEGG<-merge(outKEGG,KEGGterms[,c(3:5,7)],by=c("ID","Description"),all=T)
outGO<-merge(outGO,GOterms[,c(2:4,6)],by=c("ID","Description"),all=T)

# IX Kim et al. GWAS -----
gwas<-read.delim("cfRNA.Kim.et.al.2024.GWAS.overlaps.txt",header = F)
gwas<-gwas[gwas$V6!=-1,]
df<-cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$Gene%in%gwas$V4,]
# Perform hypergeometric test to evaluate significance of overlap
k<-195
n<-as.numeric(length(df$hgnc_symbol))
phyper(n-1, k, K-k, N, lower.tail = F)

# X GWAS browser -----
gb<-read.csv("locusbrowser.csv")
gb<-gb[gb$Gene%in%cfrna$hgnc_symbol[cfrna$Gene%in%gwas$V4],c(6,7)]
gb<-aggregate(Conclusion ~ Gene, data = gb, max)

# XI Guen et al. XWAS -----
xwas<-read.delim("cfRNA.Guen.et.al.XWAS.overlap.txt",header = F)
xwas<-xwas[xwas$V6!=-1,]
df<-cfrna[!is.na(cfrna$MetaPadj)&cfrna$MetaPadj<0.05&cfrna$Gene%in%xwas$V4,]
# Perform hypergeometric test to evaluate significance of overlap
k<-7
n<-as.numeric(length(df$hgnc_symbol))
phyper(n-1, k, K-k, N, lower.tail = F)
