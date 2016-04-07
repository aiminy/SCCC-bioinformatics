#!/usr/bin/env Rscript
#Usage: Rscript Run_DESeq_3.R inputfile_count inputfile_sample_information outputfile_prefix

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

Function.Get.DE.FC<-function(i){
count.4.sample.control<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c("control")),1]
count.4.sample.treatment<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c("treatment")),1]
index.control=which(as.character(colnames(data.count)) %in% c(as.character(count.4.sample.control)))
index.treatment=which(as.character(colnames(data.count)) %in% c(as.character(count.4.sample.treatment)))
data.count.reformated<-cbind(data.count[,index.control],data.count[,index.treatment])
condition <- factor(c("control","control","treatment","treatment"))
dds <- DESeqDataSetFromMatrix(data.count.reformated, DataFrame(condition), ~ condition)
re<-results(DESeq(dds))
re.FC<-cbind(as.data.frame(re),2^re[,2])
colnames(re.FC)[7]="FoldChange"
write.csv(cbind(re.FC,as.data.frame(counts(dds))),file=paste(output.file,"_",data.sample[i],"_DE.csv",sep=""))
write.csv(as.data.frame(colData(dds)),file=paste(output.file,"_",data.sample[i],"_info.csv",sep=""))
}

lapply(sample.index, Function.Get.DE.FC)
