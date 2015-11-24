#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
if (length(args)<=1){
	stop("Missing arguments\ninput directory, project_name, output directory, index file name, gene annotation\n**input pair-end, .txt file, with file name ends with 1 or 2 indicating the first or second of the pair**\n", call.=F)
}else if (length(args)==2){
	stop("No default output directory\n", call.=F)
}else if (length(args)==3){
	stop("No index files\n", call.=F)
}else if (length(args)==4){
	warning("No gene annotation supplied\n")
	ga = ""
}else{
	ga = args[5]
}
log_file = "run_tophat_log.txt"
cat("log output to ", log_file, "\n")

# parameters
Ncores = 4

# cread a table of inputs and run tophat alignment
# find command name in the folder
file_dir = args[1]
proj_name = args[2]
out_dir = args[3]
files_names = dir(file_dir)
files_names = grep(proj_name, files_names, value=T) 
if (length(files_names)%%2!=0){
	warning("sequence files missing at least in one pair\n")
}

file_names_unique = unique(gsub("[1-2].txt","",files_names))
P1 = paste0(file_dir, file_names_unique, "1.txt")
P2 = paste0(file_dir, file_names_unique, "2.txt")
out_dir1 = paste0(out_dir, "/", file_names_unique, "_1")
out_dir2 = paste0(out_dir, "/", file_names_unique, "_2")
table_cmd = data.frame(command = paste("tophat2 --library-type fr-unstranded -p", Ncores, "-G"), Annot=ga, 
	Outfile=paste("-o", out_dir), BowtieIndex = args[4], P1, P2)


write.table(table_cmd, "table_tophat.txt", sep=" ", quote=F, row.names=F, col.names=F)


#for(i in 1:length(P1)){
#	cmd = paste(table_cmd[i,], collapse=" ")
#	system(paste("mkdir", out_dir1[i]))
#	system(paste("mkdir", out_dir2[i]))
#	if (i==1){
#		cat("Processing ", file_names_unique[i], "\n", file=log_file)
#	}else{
#		cat("Processing ", file_names_unique[i], "\n", file=log_file, append=T)
#	}
#	cat(system(cmd, wait=F, intern=T), file=log_file, append=T)
#}


