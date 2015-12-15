**Pipeline for RNA seq data analysis**

1. Sequence short read file processing
2. Sequence alignment
3. Feature counting
   * Differential expression
  
```Rscript 
Rscript Run_DESeq_3.R inputfile_count inputfile_sample_information outputfile_prefix
```
   * Gene Set Enrichment Analysis(GSEA)
