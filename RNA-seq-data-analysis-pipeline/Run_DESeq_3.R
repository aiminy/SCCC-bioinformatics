#!/usr/bin/env Rscript
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

args = commandArgs(trailingOnly=T)

input.file = args[1]
sample.info.file=args[2]
output.file = args[3]

cat(input.file,"\n")
cat(sample.info.file,"\n")
data.count<-read.table(input.file,sep="\t",header = T,row.names=1)
data.info<-read.table(sample.info.file,sep="\t",header = F)

data.sample<-unique(trimws(as.character(data.info[,2])))
data.number.sample<-length(data.sample)

Function.Check.condition<-function(i){
sample.condition<-unique(trimws(as.character(data.info[which(data.info[,2] %in% data.sample[i]),3])))
return(sample.condition)
}

sample.index<-as.list(c(seq(1,data.number.sample)))
print(sample.index)
cat(dim(data.count),"\n")
cat(colnames(data.count),"\n")

#print(head(data.count))
#print(data.info)
#print(data.sample)
#print(data.number.sample)
#print(sample.condition)
#print(lapply(sample.index,Function.Check.condition))

#print()

Function.Get.DE.FC<-function(i){
#for (i in 1:3){
count.4.sample.control<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c("control")),1]
count.4.sample.treatment<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c("treatment")),1]

#print(as.character(count.4.sample.control))
#print(as.character(count.4.sample.treatment))
#print(which(as.character(colnames(data.count)) %in% c(as.character(count.4.sample))))
index.control=which(as.character(colnames(data.count)) %in% c(as.character(count.4.sample.control)))
index.treatment=which(as.character(colnames(data.count)) %in% c(as.character(count.4.sample.treatment)))
data.count.reformated<-cbind(data.count[,index.control],data.count[,index.treatment])
#print(head(data.count.reformated))
#data.info[which(data.info[,2] %in% data.sample[1]),3]
condition <- factor(c("control","control","treatment","treatment"))
#dds <- DESeqDataSetFromMatrix(countData.cbpK0.2)
dds <- DESeqDataSetFromMatrix(data.count.reformated, DataFrame(condition), ~ condition)
#print(colData(dds))
#dds$condition
#dds$condition <- factor(dds$condition, levels=c(rep("untreated",2),rep("treated",2)))
re<-results(DESeq(dds))
re.FC<-cbind(as.data.frame(re),2^re[,2])
#print(colnames(re.FC))
colnames(re.FC)[7]="FoldChange"
#print(colnames(re.FC))

#print(head(data.count.reformated))
#print(counts(dds))
#print(names(re))
#re.count<-merge(re,data.count.r,by=0)
#re.treatment.control=merge(re.treatment,count.4.sample.control,by=0)
write.csv(cbind(re.FC,as.data.frame(counts(dds))),file=paste(output.file,"_",data.sample[i],"_DE.csv",sep=""))
write.csv(as.data.frame(colData(dds)),file=paste(output.file,"_",data.sample[i],"_info.csv",sep=""))
}

lapply(sample.index, Function.Get.DE.FC)