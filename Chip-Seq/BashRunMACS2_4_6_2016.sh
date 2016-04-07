#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J macs2_chip_seq
#BSUB -P macs2_chip_seq
#BSUB -o %J.out
#BSUB -e %J.err

mkdir /scratch/projects/bbc/Peak_chip_seq

while read line; do

#f1=`echo "$line" | awk -F"\t" '{print $1}'`                                                                                                                                                                
#f2=`echo "$line" | awk -F"\t" '{print $2}'`                                                                                                                                                                

f=`echo "$line"`

echo "$f"

dir_name=$(dirname "$f")
file_name=$(basename "$f")

echo "$dir_name"
echo "$file_name"

sample_name=`echo "$file_name" | awk -F"." '{print $1}'`

echo "$sample_name"

cat > ~/Script_bash/Run_"$sample_name"_to_chip_seq_peak.sh <<EOF

#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J macs2_chip_seq
#BSUB -P macs2_chip_seq
#BSUB -o %J.out
#BSUB -e %J.err

python /nethome/axy148/MACS/bin/macs2 callpeak -t /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position.bam -f BAM -g mm -n "$sample_name"_bam_mm --outdir /scratch/projects/bbc/Peak_chip_seq -B -q 0.01

EOF

bsub -e %J.err -P Bioinformatics4count < ~/Script_bash/Run_"$sample_name"_to_chip_seq_peak.sh

done < $1
