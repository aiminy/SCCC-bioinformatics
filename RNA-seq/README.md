**Pipeline for RNA seq data analysis**

* Download this repository
```
git clone https://github.com/aiminy/SCCC-bioinformatics.git
```

* Short read sequence files processing
```
#download FastQC, then type
mkdir QC_output
fastqc -o QC_output path/*fastq.gz

#move to the directory(QC_output) containg the files ending with _fastqc.zip generated with FastQC
#then run: 
python ~/SCCC-bioinformatics/RNA-seq-data-analysis-pipeline/QcSum.py QC_summary_output_file
```
  * Sequence alignment
```bash
#To process for pair end files,run the following:
sh ~/SCCC-bioinformatics/RNA-seq-data-analysis-pipeline/RunAlignment4PairEnd.sh Sample_9.txt
#To process for single end files, run the following:
sh ~/SCCC-bioinformatics/RNA-seq-data-analysis-pipeline/RunAlignment4SingleEnd.sh Sample_9.txt  
```
   * Feature counting
```bash
# put all bam file names into FileBam.txt
sh ~/SCCC-bioinformatics/RNA-seq-data-analysis-pipeline/Bash_run_submit_job_4_sorted_bam.sh FileBam.txt

#put all sorted bam names into FileBamSorted.txt 
sh ~/SCCC-bioinformatics/RNA-seq-data-analysis-pipeline/Bash_run_submit_job_4_sorted_bam_2.sh FileBamSorted.txt
```

```python
#To put all count files into one file
python ~/SCCC-bioinformatics/RNA-seq-data-analysis-pipeline/merge_tables.py koji_data_count_sample_information_c1_c3_3.txt Output_count_samples24.txt
```

   * Differential expression
  
```Rscript 
Rscript ~/SCCC-bioinformatics/RNA-seq-data-analysis-pipeline/Run_DESeq_3.R inputfile_count inputfile_sample_information outputfile_prefix

```
   * Gene Set Enrichment Analysis(GSEA)
