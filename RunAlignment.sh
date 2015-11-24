#! /bin/bash

#paste <(ls *1.txt* | sort  -t"_" -k3)  <(ls *2.txt* | sort  -t"_" -k3) > Fq_24.txt

#<<EOF
while read line; do

f1=`echo "$line" | awk -F"\t" '{print $1}'`
f2=`echo "$line" | awk -F"\t" '{print $2}'`

#f=`echo "$line"`

echo "$f1\t$f2"
sample_name=`echo "$f1" | awk -F"_" '{print $3}'`

#echo "$sample_name"_aligned_reads.sam

'bowtie -p 8 --chunkmbs 2000 --sam mm10_index/genome -1 "$f1" -2 "$f2" "$sample_name"_aligned_reads.sam

done < $1
#EOF

bowtie -p 8 --chunkmbs 2000 --sam mm10_index/genome -1 <(gunzip nBishopric_Project1_201348190-01_S_1_1.txt.gz) -2 <(gunzip nBishopric_Project1_201348190-01_S_1_2.txt.gz)  201348190-01_aligned_reads.sam&




