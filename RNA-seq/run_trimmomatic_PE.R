#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
if (length(args)<1){
	stop("Usage: Rscript run_trimmomatic_SE input directory, output directory, adapter file (.fa) \n**input single-end file**\n", call.=F)
}else if (length(args)==1){
	stop("No default output directory\n", call.=F)
}else if (length(args)==2){
	warning("No adapter seq\n", call.=F)
	adpClip = ""
}

# parameters
Ncores = 8

# cread a table of inputs
# find command name in the folder
file_dir = args[1]
out_dir = args[2]
file_names = dir(file_dir)

if(!exists("adpClip")){
	adpClip = paste0("ILLUMINACLIP:",args[3],":2:30:10")
}
leading = "LEADING:3"
trailing = "TRAILING:3"
slindingwindow = "SLIDINGWINDOW:4:20"
minlen = "MINLEN:20"
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
out_dir1 = paste0(out_dir, "/trimmed_", file_names_unique, "1", ext_pttn)
out_dir2 = paste0(out_dir, "/trimmed_", file_names_unique, "2", ext_pttn)

table_cmd = cbind(command = paste("java -jar /share/apps/trimmomatic/0.32/trimmomatic-0.32.jar PE -threads", Ncores, "-phred33"), P1, P2,Outfile=paste0(out_dir1," ", out_dir1, ".unpaired ", out_dir2, " ", out_dir2, ".unpaired"), op1 = adpClip, op2 = leading, op3 = trailing, op4 = slindingwindow, op5 = minlen)

#write.table(table_cmd, "table_trimmer.txt", sep=" ", quote=F, row.names=F, col.names=F)

out_bsubscript = "bsub_trimPE_tmp.txt"
cat(
paste("#!/bin/bash\n#BSUB -J bsub_trim", "\n#BSUB -o %J.trimming.log\n#BSUB -e %J.trimming.err\n#BSUB -W 6:00\n#BSUB -q general\n#BSUB -n", Ncores, "\n#BSUB -R 'rusage[mem=4096] span[hosts=1]'\n#BSUB -N\n#BSUB -u yban@med.miami.edu"), file=out_bsubscript)
cat("module load java/1.8.0_60\n", file=out_bsubscript, append=T)
cat("module load trimmomatic/0.32\n", file=out_bsubscript, append=T)


for(i in 1:length(P1)){
	cat(paste(table_cmd[i,], collapse=" "),"\n", file=out_bsubscript, append=T)
}

system(paste("bsub <",out_bsubscript))
system(paste("rm", out_bsubscript))

