library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library(biomaRt)

out_gsea = "/media/MyDATA/Result/Postanalysis/GSEA/Bishopric/"

#########################################################################################
# Description (main)
# function takes a merged count table with one comparison,  and pvalues, foldchanges
# output: if less than 2 samples use prerank for GSEA, if 3 or more use VST normalization and outpur matrix to run GSEA
############################################################################################


# input a merged count table


#----- clean data -----#
# read csv input count table
#  and then normaliz, output to GSEA format
# remove summary data from HTSeq count table
removeSF = function(data_table){
  Fname = c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
  Fname = c(Fname, toupper(Fname))
  Temp = na.omit(match(Fname, rownames(data_table)))
  if(length(Temp)==0){
    return(data_table)
  }
  return(data_table[-Temp, ])
}


# define experiment comparisons


# run DESeq analysis
# DESeq_2group: comparing trt with control
# DESeq_mlutigroupLRT: overall test whether there is difference among groups
DESeq_2group = function(Data, Sample_names1, Sample_names2){
  cat(paste("Comparing", Sample_names1, "(Treated)\n\tto", Sample_names2, "(Untreated)\n"))
  N_data = nrow(Data)
  nsample1 = length(Sample_names1)
  nsample2 = length(Sample_names2)
  colID_select = match(c(Sample_names1,Sample_names2), colnames(Data))
  Data_s = Data[, colID_select]
  condition = factor(c(rep("Treated",nsample1),rep("Untreated",nsample2)))
  dds <- DESeqDataSetFromMatrix(Data_s, DataFrame(condition), ~ condition)
  return(results(DESeq(dds)))
}

# Example: test.de = DESeq_2group(data.cbpKO, colnames(data.cbpKO)[1:2], colnames(data.cbpKO)[3:4])

DESeq_mlutigroupLRT = function(SampleINFO, count_DIR){
  sample.info = read.table(paste0(count_DIR,"/",SampleINFO),sep="\t",as.is=T)
  colnames(SampleINFO) = c("sampleName","fileName","condition")
  dds1 = DESeqDataSetFromHTSeqCount(sampleTable=SampleINFO,directory=count_DIR,design= ~ condition)
  data.count.3.groups<-counts(dds1)
  ddsLRT <- DESeq(dds1, test="LRT", reduced= ~ 1)
  return(results(ddsLRT))
}
# Example: DESeq_mlutigroupLRT("Sample_9_2.txt", "/media/MyDATA/Result/Count/Bishopric/")


# organize output result pvalue and count table


# generate GSEA inputs
# output: if less than 2 samples use prerank for GSEA, if 3 or more use VST normalization and outpur matrix to run GSEA
#----for more than 2 samples, use function GSEA_twoclass----#
GSEA_twoclass=function(Data, Sample_names1, Sample_names2, outfile_name, out_gsea){
  N_data = nrow(Data)
  nsample1 = length(Sample_names1)
  nsample2 = length(Sample_names2)
  colID_select = match(c(Sample_names1,Sample_names2), colnames(Data))
  Data_s = Data[, colID_select]
  Data_s<-cbind(Name=rownames(Data),description="NA", Data_s)
  name1=outfile_name
  name2="gct"
  filename=paste(out_gsea,name1,".",name2,sep="")
  cat("#1.2\n",file=filename,append=FALSE)
  secondline=paste(as.character(N_data),as.character(nsample1+nsample2),sep="\t")
  cat(secondline,file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  write.table(Data_s,file=filename,append=TRUE,row.names=FALSE,col.name=TRUE,sep="\t", quote=FALSE)
  
  class.bin = t(c(rep(0,nsample1),rep(1,nsample2)))
  name3="cls"
  filename1=paste(out_gsea,name1,".",name3,sep="")
  cat(paste(nsample1+nsample2,"2 1\n"),file=filename1,append=FALSE)
  secondline1=paste("# high low","\n")
  cat(secondline1,file=filename1,append=TRUE)
  write.table(class.bin,file=filename1,append=TRUE,row.names=FALSE,col.name=FALSE)
}

# use example : 
# GSEA_twoclass(vsd.cbpkO, colnames(vsd.cbpkO)[1:2], colnames(vsd.cbpkO)[3:4], "cbpKO", out_gsea1)

#----------------------------------#
# use preranked for GSEA
#----------------------------------#

# read csv input DESeq result table
#  and then output to GSEA preranked format
diff.cbpKO = read.csv("/media/MyDATA/Result/MergedCount/Bishopric/DE_count_cbpK0.csv", row.names=1)
rownames(diff.cbpKO) = toupper(rownames(diff.cbpKO))
diff.p300K0 = read.csv("/media/MyDATA/Result/MergedCount/Bishopric/DE_count_p300K0.csv", row.names=1)
rownames(diff.p300K0) = toupper(rownames(diff.p300K0))
diff.WT = read.csv("/media/MyDATA/Result/MergedCount/Bishopric/DE_count_WT.csv", row.names=1)
rownames(diff.WT) = toupper(rownames(diff.WT))

diff.cbpKO = removeSF(diff.cbpKO)
diff.p300K0 = removeSF(diff.p300K0)
diff.WT = removeSF(diff.WT)

rnk.cbpKO = na.omit(cbind(rownames(diff.cbpKO), -diff.cbpKO$log2FoldChange))
rnk.p300K0 = na.omit(cbind(rownames(diff.p300K0), -diff.p300K0$log2FoldChange))
rnk.WT = na.omit(cbind(rownames(diff.WT), -diff.WT$log2FoldChange))

write.table(rnk.cbpKO, paste0(out_gsea1, "cbpKO_preranked.rnk"),row.names=F,col.names=F,
            sep="\t",quote=F)
write.table(rnk.p300K0, paste0(out_gsea1, "p300K0_preranked.rnk"),row.names=F,col.names=F,
            sep="\t",quote=F)
write.table(rnk.WT, paste0(out_gsea1, "WT_preranked.rnk"),row.names=F,col.names=F,
            sep="\t",quote=F)