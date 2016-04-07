#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
home_mm10_index = "/nethome/yxb173/Genome_Ref/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
home_mm10_ga = "/nethome/yxb173/Genome_Ref/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
if (length(args)<1){
	stop("Usage(Pair-end RNASeq): Rscript alignment2count_PE.R input fastq directory, output alignment directory, output count directory, index file name, gene annotation\n", call.=F)
}else if (length(args)==1){
	stop("No default output alignment directory\n", call.=F)
}else if (length(args)==2){
	stop("No default output count directory\n", call.=F)
}else if (length(args)==3){
	warning(paste("No index files, set default:\n",home_mm10_index), call.=F)
}else if (length(args)==4){
	warning(paste("No gene annotation supplied, set default:\n",home_mm10_ga), call.=F)
}

# parameters
Ncores = 8

# cread a table of inputs and run tophat alignment
# find command name in the folder
file_dir = args[1]
align_dir = args[2]
count_dir = args[3]
if(is.na(args[4])){
	indexf = home_mm10_index
}else{
	indexf = args[4]
}
if(is.na(args[5])){
	ga = home_mm10_ga
}else{
	ga = args[5]
}

file_names = dir(file_dir)
if (length(file_names)%%2!=0){
	warning("sequence files missing at least in one pair\n")
}

ext_pttn = c(".txt", ".fastq")
tmp = NULL
for (p in ext_pttn){tmp = c(tmp, any(grepl(p,file_names)))}
if(!any(tmp)){
	stop("Cannot detect file extension format, please modifiy ext_pttn\n")
}else if(length(ext_pttn[tmp])>1){
	warning("Multiple file extension formats detected, using the first one\n")
}

ext_pttn = ext_pttn[tmp][1]

file_names_unique = unique(gsub(paste0("[1-2]", ext_pttn),"",file_names))
P1 = paste0(file_dir, "/", file_names_unique, "1", ext_pttn)
P2 = paste0(file_dir, "/", file_names_unique, "2", ext_pttn)

out_dir1 = paste0(align_dir, "/", file_names_unique)

table_cmd = cbind(command = paste("tophat2 --library-type fr-unstranded -p", Ncores, "-G"), Annot=ga, Outfile=paste("-o", out_dir1), BowtieIndex = indexf, P1, P2)

table_cmd2 = cbind(com="samtools sort -n", file=paste0(out_dir1,"/accepted_hits.bam"), out=paste0(out_dir1,"/accepted_hits.sorted"))

table_cmd3 = cbind(com="htseq-count -f bam -r name -s no", file=paste0(out_dir1,"/accepted_hits.sorted.bam"), gtf_dir=ga, out=paste0("> ",count_dir,"/",file_names_unique,".count"))

#write.table(table_cmd, "table_tophat.txt", sep=" ", quote=F, row.names=F, col.names=F)

for(i in 1:length(P1)){ 
	out_bsubscript = paste0(file_names_unique[i],".txt")
	cat(
paste("#!/bin/bash\n#BSUB -J", file_names_unique[i], "\n#BSUB -o %J.align2count.log\n#BSUB -e %J.align2count.err\n#BSUB -W 6:00\n#BSUB -q general\n#BSUB -n 8\n#BSUB -R 'rusage[mem=4096] span[hosts=1]'\n#BSUB -N\n#BSUB -u\n"), file=out_bsubscript)
	cat("module load bowtie2\n", file=out_bsubscript, append=T)
	cat("module load tophat/2.0.11\n", file=out_bsubscript, append=T)
	cat("module load samtools/0.1.19\n", file=out_bsubscript, append=T)
	cat(paste(table_cmd[i,], collapse=" "),"\n", file=out_bsubscript, append=T)
        cat(paste(table_cmd2[i,], collapse=" "),"\n", file=out_bsubscript, append=T)
        cat(paste(table_cmd3[i,], collapse=" "), "\n", file=out_bsubscript, append=T)

	system(paste("mkdir", out_dir1[i]))
	system(paste("bsub <",out_bsubscript))
	system(paste("rm", out_bsubscript))
}


