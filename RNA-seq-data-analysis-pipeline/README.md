**Pipeline for RNA seq data analysis**

*Short read sequence files processing
```
#download FastQC, then type
mkdir QC_output
fastqc -o QC_output fq1 fq2 fq3

#move to the directory(QC_output) containg the files ending with _fastqc.zip generated with FastQC
#then run: 
python ~/Code/RNA-seq-data-analysis-pipeline/QcSum.py QC_summary_output_file
```
  * Sequence alignment
```bash
#To process for pair end files,run the following:
sh ~/Script_bash/RunAlignment4PairEnd.sh Sample_9.txt
#To process for single end files, run the following:
sh ~/Script_bash/RunAlignment4SingleEnd.sh Sample_9.txt  
```
   * Feature counting
```bash
# put all bam file names into FileBam.txt
sh Bash_run_submit_job_4_sorted_bam.sh FileBam.txt

#put all sorted bam names into FileBamSorted.txt 
sh Bash_run_submit_job_4_sorted_bam_2.sh FileBamSorted.txt
```

```python
#To put all count files into one file
python python ~/Code/RNA-seq-data-analysis-pipeline/merge_tables.py Sample_10_raw_count.txt Output_count_samples10.txt

```

   * Differential expression
  
```Rscript 
Rscript Run_DESeq_3.R inputfile_count inputfile_sample_information outputfile_prefix
```
   * Gene Set Enrichment Analysis(GSEA)
