library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")
library(gplots)
library(RColorBrewer)
library(survival)
library(limma)
library(edgeR)
library(multtest)
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
biocLite("goseq")
library(goseq)

setwd("~/bishopric/RNAseq/rawData/140128_SN392_0190_AC3EYUACXX")
data.140128.sn392.0190<-read.table("merged_counts.txt",sep="\t",header = T,row.names=1)

#For genotype cbpKO
data.cbpKO<-data.140128.sn392.0190[,c(1,9,2,3)]
head(data.cbpKO)

countData <- data.cbpKO
condition <- factor(c("cbpKO_t","cbpKO_t","cbpKO_ck","cbpKO_ck"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
re.cbpKO<-DESeq(dds)
str(re.cbpKO)
save.image(file="DE_cbpKO.RData")
dir(pattern = "RData")















