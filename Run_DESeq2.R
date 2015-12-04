source("http://bioconductor.org/biocLite.R")
biocLbiocLite("goseq")
bioite("biomaRt")

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
library(biomaRt)
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
write.csv(data.WT,file="Count_WT.csv")

#For genotype p300K0
data.p300K0<-data.140128.sn392.0190[,c(8,9,10,11)]
head(data.p300K0)
countData.p300K0 <- data.p300K0
condition.p300K0 <- factor(c("p300K0_t","p300K0_t","p300K0_ck","p300K0_ck"))
dds.p300K0 <- DESeqDataSetFromMatrix(countData.p300K0, DataFrame(condition.p300K0), ~ condition.p300K0)
re.p300K0<-results(DESeq(dds.p300K0))
write.csv(re.p300K0,file="DE_p300K0.csv")
write.csv(data.p300K0,file="Count_p300K0.csv")

#For Jianping data
data.jianping<-read.table("merged_counts.txt",sep="\t",header = T,row.names=1)
head(data.jianping)

##For HET vs K
data.het.k<-data.jianping[,1:6]
head(data.het.k)
countData <- data.het.k
condition <- factor(c("HET","HET","HET","K","K","K"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
re.het.k<-results(DESeq(dds))
write.csv(re.het.k,file="DE_HET_K_2.csv")
write.csv(data.het.k,file="Count_HET_K_2.csv")
write.csv(cbind(re.het.k,data.het.k),file="DE_count_HET_K_2.csv")

#For HET vs WT
data.het.wt<-data.jianping[,c(1:3,7:9)]
head(data.het.wt)
countData <- data.het.wt
condition <- factor(c("HET","HET","HET","WT","WT","WT"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
re.het.wt<-results(DESeq(dds))
write.csv(re.het.wt,file="DE_HET_WT_2.csv")
write.csv(data.het.wt,file="Count_HET_WT_2.csv")
write.csv(cbind(re.het.k,data.het.k),file="DE_count_HET_WT_2.csv")

#For K vs WT
data.k.wt<-data.jianping[,4:9]
head(data.k.wt)
countData <- data.k.wt
condition <- factor(c("K","K","K","WT","WT","WT"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
re.k.wt<-results(DESeq(dds))

head(re.k.wt)
head(data.k.wt)

write.csv(re.k.wt,file="DE_K_WT_2.csv")
write.csv(data.k.wt,file="Count_K_WT_2.csv")
write.csv(cbind(re.k.wt,data.k.wt),file="DE_count_K_WT_2.csv")



samples <- read.table("Sample_9_2.txt",sep="\t",as.is=T)
samples
colnames(samples)<-c("sampleName","fileName","condition")
dds1 <- DESeqDataSetFromHTSeqCount(sampleTable=samples,directory="./",design= ~ condition)
data.count.3.groups<-counts(dds1)
class(data.count.3.groups)

which(!rownames(data.jianping) %in% rownames(data.count.3.groups))
data.jianping[c(4470,5189,5807,10727,12906),]
dim(data.count.3.groups)

dds1 <- DESeq(dds1) 
res <- results(dds1)
res

res.HET.K <- results(dds1, contrast=c("condition","HET","K"))
res.K.HET <- results(dds1, contrast=c("condition","K","HET"))
res.HET.K
res.K.HET

res.HET.WT <- results(dds1, contrast=c("condition","HET","WT"))
res.WT.HET <- results(dds1, contrast=c("condition","WT","HET"))
res.HET.WT
res.WT.HET

res.K.WT <- results(dds1, contrast=c("condition","K","WT"))
res.WT.K <- results(dds1, contrast=c("condition","WT","K"))
res.K.WT
res.WT.K

ddsLRT <- DESeq(dds1, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)


#dim(resLRT)
#class(resLRT)
#names(resLRT)
dim(data.jianping)
dim(resLRT)
dim(data.count.3.groups)

resLRT.data.count.3.groups<-cbind(resLRT,data.count.3.groups)
resLRT.data.jianping<-cbind(resLRT,data.jianping[-c(4470,5189,5807,10727,12906),])
head(resLRT)
head(data.jianping[-c(4470,5189,5807,10727,12906),])

which(rownames(data.jianping) %in% c("0610005C13Rik"))

data.jianping[8960,]

resLRT.data.jianping.2<-merge(as.data.frame(resLRT),data.jianping[-c(4470,5189,5807,10727,12906),],by=0)

head(resLRT.data.jianping.2)
head(resLRT.data.jianping)
dim(resLRT.data.jianping.2)


st(resLRT)
write.csv(cbind(res.HET.K,res.HET.WT,res.K.WT,resLRT),file="All_plus_contrast_comparision.csv")
write.csv(resLRT,file="All_3_groups_comparision.csv")
write.csv(resLRT.data.jianping.2,file="All_3_groups_raw_count.csv")

#str(re.cbpKO)
#y<-DGEList(counts=data.cbpKO, genes=rownames(data.cbpKO))
#names(y)
#head(y$counts)
#head(y$samples)
#tail(y$genes)
#y<-calcNormFactors(y)

#For genotype WT
#data.WT<-data.140128.sn392.0190[,c(4,5,6,7)]
#head(data.WT)
#countData.WT <- data.WT
#condition.WT <- factor(c("WT_t","WT_t","WT_ck","WT_ck"))
#dds.WT <- DESeqDataSetFromMatrix(countData.WT, DataFrame(condition.WT), ~ condition.WT)
#re.WT<-results(DESeq(dds.WT))
#results(re.WT)
#str(re.WT)
#write.csv(re.WT,file="DE_WT.csv")
#write.csv(data.cbpKO,file="Count_WT.csv")

#For genotype p300K0
#data.p300K0<-data.140128.sn392.0190[,c(8,9,10,11)]
#head(data.p300K0)
#countData.p300K0 <- data.p300K0
#condition.p300K0 <- factor(c("p300K0_t","p300K0_t","p300K0_ck","p300K0_ck"))
#dds.p300K0 <- DESeqDataSetFromMatrix(countData.p300K0, DataFrame(condition.p300K0), ~ condition.p300K0)
#re.p300K0<-results(DESeq(dds.p300K0))
#write.csv(re.WT,file="DE_p300K0.csv")
#write.csv(data.cbpKO,file="Count_p300K0.csv")

#str(re.p300K0)
#rld<-rlog(dds.p300K0)
#plotPCA(rld)
#counts(dds.p300K0)
#assay(dds.p300K0)
#fpm(dds.p300K0)
#fpkm(dds.p300K0)
#re.p300K0<-DESeq(dds.p300K0)
#results(re.p300K0)
#str(re.p300K0)

#save.image(file="DE_cbpKO.RData")
#dir(pattern = "RData")