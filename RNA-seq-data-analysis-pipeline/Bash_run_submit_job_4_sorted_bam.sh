#! /bin/bash

#sh Bash_run_submit_job_4_sorted_bam.sh FileBam.txt
while read line; do

f=`echo "$line"`

echo "$f"

sample_name=`echo "$f" | awk -F"." '{print $1}'`

echo "$sample_name"

cat > Run_"$sample_name"_htseq_count_4_sorted_bam.sh <<EOF

samtools sort -n "$f" "$f"_sorted.bam

#htseq-count -f bam -r name -s yes -i gene_name "$f"_sorted.bam genes.gtf  > "$sample_name"_raw_count_4_sorted_bam.txt

EOF

bsub -P Bioinformatics4count < Run_"$sample_name"_htseq_count_4_sorted_bam.sh

done < $1