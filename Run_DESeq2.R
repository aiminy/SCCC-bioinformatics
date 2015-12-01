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
data.cbpKO<-data.140128.sn392.0190[,c(1,12,2,3)]
head(data.cbpKO)

countData <- data.cbpKO
condition <- factor(c("cbpKO_t","cbpKO_t","cbpKO_ck","cbpKO_ck"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
re.cbpKO<-results(DESeq(dds))
write.csv(re.cbpKO,file="DE_cbpK0.csv")
write.csv(data.cbpKO,file="Count_cbpk0.csv")

#str(re.cbpKO)
#y<-DGEList(counts=data.cbpKO, genes=rownames(data.cbpKO))
#names(y)
#head(y$counts)
#head(y$samples)
#tail(y$genes)
#y<-calcNormFactors(y)

#For genotype WT
data.WT<-data.140128.sn392.0190[,c(4,5,6,7)]
head(data.WT)
countData.WT <- data.WT
condition.WT <- factor(c("WT_t","WT_t","WT_ck","WT_ck"))
dds.WT <- DESeqDataSetFromMatrix(countData.WT, DataFrame(condition.WT), ~ condition.WT)
re.WT<-results(DESeq(dds.WT))
#results(re.WT)
#str(re.WT)
write.csv(re.WT,file="DE_WT.csv")
write.csv(data.cbpKO,file="Count_WT.csv")

#For genotype p300K0
data.p300K0<-data.140128.sn392.0190[,c(8,9,10,11)]
head(data.p300K0)
countData.p300K0 <- data.p300K0
condition.p300K0 <- factor(c("p300K0_t","p300K0_t","p300K0_ck","p300K0_ck"))
dds.p300K0 <- DESeqDataSetFromMatrix(countData.p300K0, DataFrame(condition.p300K0), ~ condition.p300K0)
re.p300K0<-results(DESeq(dds.p300K0))
write.csv(re.WT,file="DE_p300K0.csv")
write.csv(data.cbpKO,file="Count_p300K0.csv")


#str(re.p300K0)




rld<-rlog(dds.p300K0)
plotPCA(rld)
counts(dds.p300K0)
assay(dds.p300K0)
fpm(dds.p300K0)
fpkm(dds.p300K0)
re.p300K0<-DESeq(dds.p300K0)
results(re.p300K0)
str(re.p300K0)

save.image(file="DE_cbpKO.RData")
dir(pattern = "RData")















