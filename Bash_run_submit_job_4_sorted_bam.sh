#! /bin/bash
#BSUB -n 32                                                                                                                                                                                                
#BSUB -q general                                                                                                                                                                                           
#BSUB -W 05:00                                                                                                                                                                                             
#BSUB -J Sort_BAM_bash                                                                                                                                                                                     
#BSUB -P Sort_Bam_328_all                                                                                                                                                                                  
#BSUB -o %J.out                                                                                                                                                                                            
#BSUB -e %J.err

#sh Bash_run_submit_job_4_sorted_bam.sh FileBam.txt
while read line; do

f=`echo "$line"`

echo "$f"

sample_name=`echo "$f" | awk -F"." '{print $1}'`

sample_name2=$(basename "$sample_name")

echo "$sample_name"

cat > ~/Script_bash/Run_"$sample_name2"_htseq_count_4_sorted_bam.sh <<EOF
#!/bin/bash                                                                     
#BSUB -n 32                                                                     
#BSUB -q general                                                                
#BSUB -W 05:00                                                                  
#BSUB -J Sort_BAM_each                                                                 
#BSUB -P Sort_BAM_each_p                                                        
#BSUB -o %J.out                                                                 
#BSUB -e %J.err     

samtools sort -n "$f" "$f"_sorted.bam

#htseq-count -f bam -r name -s yes -i gene_name "$f"_sorted.bam genes.gtf  > "$sample_name"_raw_count_4_sorted_bam.txt

EOF

bsub -P Bioinformatics4count < ~/Script_bash/Run_"$sample_name2"_htseq_count_4_sorted_bam.sh

done < $1