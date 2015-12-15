**Pipeline for RNA seq data analysis**

1. Sequence short read file processing
2. Sequence alignment

```bash
#To process for pair end files,run the following:
sh ~/Script_bash/RunAlignment4PairEnd.sh Sample_9.txt
#To process for single end files, run the following:
sh ~/Script_bash/RunAlignment4SingleEnd.sh Sample_9.txt  
```

3. Feature counting
   * Get counts for each gene
```bash
# put all bam file names into FileBam.txt
sh Bash_run_submit_job_4_sorted_bam.sh FileBam.txt

#put all sorted bam names into FileBamSorted.txt 
sh Bash_run_submit_job_4_sorted_bam_2.sh FileBamSorted.txt
```
   * Differential expression
  
```Rscript 
Rscript Run_DESeq_3.R inputfile_count inputfile_sample_information outputfile_prefix
```
   * Gene Set Enrichment Analysis(GSEA)
