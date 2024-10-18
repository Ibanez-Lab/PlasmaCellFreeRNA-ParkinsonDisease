#!/usr/bin/env Rscript

# 0 Set up the environment -----
## Load the necessary libraries ----
library(tximport)
library(matrixStats)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(dplyr)
library(data.table)
## Define a "not in" opperator ----
`%!in%` <- Negate(`%in%`)

# Use the following color-blind-friendly colours for plots whenever possible
myPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999","#0072B2", "#D55E00", "#CC79A7","#000000" )

# 1 Load the data -----
# Import salmon quantification files
quantFiles <- list.files("data_path/",pattern="quant.genes.sf",recursive = TRUE,full.names = TRUE)
name<-list.files("data_path")/
names(quantFiles) <- name
quantData <- tximport(quantFiles,txOut=TRUE,type="salmon")
countMatrix <- quantData$counts
#save the count matrix for subsequent analyses
write.csv(countMatrix,"cfRNA.count.matrix.csv")

# Load Picard metrics table
summaryMetrics <- read.delim("multiqc_picard_AlignmentSummaryMetrics.txt", sep = "\t")
summaryMetrics %>% select(-(SAMPLE:READ_GROUP)) 

# Load metadata for the samples
metadata <- read.csv("sample_metadata.csv")

# Merge summary stats & metadata
sumMergeDf <- merge(summaryMetrics, metadata, by.x = "Sample", "SampleID")

# Remove any genes with less than 10 counts in 90% or more of samples
countMatrix <- countMatrix[(rowCounts(countMatrix<10) < round(0.9*dim(countMatrix)[2])),] #removed 32301 genes

# Ensure that summary metrics and metadata are available for all samples
countMatrix<-countMatrix[,colnames(countMatrix)%in%rownames(sumMergeDf)]
sumMergeDf<-sumMergeDf[rownames(sumMergeDf)%in%colnames(countMatrix),]

# 2 Quality Control -----
rows<-row.names(countMatrix)
countMatrix<-sapply(countMatrix, as.integer)
row.names(countMatrix)<-rows
dds<-DESeqDataSetFromMatrix(countData = countMatrix, colData = sumMergeDf, design = ~ 1)

## Normalize and log transform the counts -----
rld <- rlog(dds, blind=TRUE)
 
## Perform PCA -----
PCA<-plotPCA(rld, intgroup = "group", returnData = T)
percentVar <- round(100 * attr(PCA, "percentVar"))
pc1.sd <- sd(pcaData$PC1)
pc1.mean <- mean(pcaData$PC1)
pc2.sd <- sd(pcaData$PC2)
pc2.mean <- mean(pcaData$PC2)

range.pc1 <- c(pc1.mean - 3*pc1.sd, pc1.mean + 3*pc1.sd)
range.pc2 <- c(pc2.mean - 3*pc2.sd, pc2.mean + 3*pc2.sd)

# check PCA plots with SD borders
panel <- c("#E69F00", "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pcaData <- pcaData %>% mutate(group = NULL)
ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size=5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_hline(yintercept = range.pc2, linetype = 'dashed') +
  geom_vline(xintercept = range.pc1, linetype = 'dashed') +
  scale_color_manual(values=panel) + 
  theme_bw()


## Correlation analyses-----
# Perform correlation analyses to see what is 
# affecting the PCA

# Combine summary metrics and PCA data frames
sumMergeDf$name<-rownames(sumMergeDf)
corDF <- merge(PCA, sumMergeDf, by = 'name')

# Select only numeric columns with no missing values
corDF <- select_if(corDF, is.numeric)
corDF <- corDF[colSums(is.na(corDF)) < 1]
corDF <- corDF[colSums(corDF) !=0]

# Perform Pearson correlation test and format 
# the correlation matrix for plotting
corMatrix <- round(cor(corDF),2)
corMatrix[lower.tri(corMatrix)] <- NA
corMatrixMelt <- melt(corMatrix, na.rm = TRUE)

# Make a heatmap showing the r values for each respective correlation
ggplot(corMatrixMelt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#E69F00", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),axis.text.y = element_text(size = 8))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank(), legend.justification = c(1, 0), legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

# Make a heatmap showing teh significance of each tested correslationdf_cor <- rcorr(as.matrix(corDF), type = "pearson")
df_long <- inner_join(melt(df_cor$r, value.name = "r"), melt(df_cor$P, value.name = "p"), by = c("Var1", "Var2")) %>%
  rowwise() %>% mutate(pair = sort(c(Var1, Var2)) %>% paste(collapse = ",")) %>% group_by(pair) %>% distinct(pair, .keep_all = T)
df_long[is.na(df_long)]<-1
df_long$sig <- "*"
df_long$sig[df_long$p>0.05]<-""
corP<-ggplot(df_long, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile() +
  coord_equal()+
  scale_fill_gradient2(low = "#0072B2", high = "#E69F00", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
  geom_text(aes(label = sig), color = "black", size = 4) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank())

## Remove outliers -----
outliers<- PCA$name[abs(PCA$PC1)>(mean(PCA$PC1)+3*sd(PCA$PC1)) | abs(PCA$PC2)>(mean(PCA$PC2)+3*sd(PCA$PC2))]
countMatrix <- countMatrix[,colnames(countMatrix) %!in% outliers]
sumMergeDf<-sumMergeDf[rownames(sumMergeDf)%in%colnames(countMatrix),]

dds<-DESeqDataSetFromMatrix(countData = countMatrix, colData = sumMergeDf, design = ~ 1)
rld <- rlog(dds, blind=TRUE)
save(rld, file = "rlog.cfRNA.counts.RData")
