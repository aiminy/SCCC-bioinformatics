#! /bin/bash

#sh ~/Script_bash/RunAlignment4PairEnd.sh Sample_9.txt  
while read line; do

f1=`echo "$line" | awk -F"\t" '{print $1}'`
f2=`echo "$line" | awk -F"\t" '{print $2}'`

echo "$f1\t$f2"
sample_name=`echo "$f1" | awk -F"_" '{print $3}'`

echo "$sample_name"

cat > Run_"$sample_name"_tophat.sh <<EOF
tophat -G genes.gtf -p 4 -o "$sample_name"_tophat_out mm10_index_bt2/genome ~/RNAseqData/"$f1" ~/RNAseqData/"$f2"
mv "$sample_name"_tophat_out/accepted_hits.bam "$sample_name".bam
EOF

done < $1
