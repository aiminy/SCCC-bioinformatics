#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J macs2_chip_seq
#BSUB -P macs2_chip_seq
#BSUB -o %J.out
#BSUB -e %J.err

#Download reference genome
#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm10.chromInfo" > mm10.genome

#! /bin/bash                                                                                                                                                                                                
#sh ~/Script_bash/BashRunConvertBamtoBigwig.sh FileBam.txt

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

#mkdir /scratch/projects/bbc/BAM2BW

cat > ~/Script_bash/Run_"$sample_name"_to_chip_seq_peak.sh <<EOF

#!/bin/bash
#BSUB -n 32
#BSUB -q general
#BSUB -W 05:00
#BSUB -J macs2_chip_seq
#BSUB -P macs2_chip_seq
#BSUB -o %J.out
#BSUB -e %J.err

#samtools sort "$f" /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position
#genomeCoverageBed -ibam /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position.bam -bg -g ~/Mus_musculus/mm10.genome > /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position.bdg
#remove header in mm10.genome, and change name to mm10_2.genome
#LC_COLLATE=C sort -k1,1 -k2,2n /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position.bdg > /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position_sorted.bdg 
#bedGraphToBigWig  /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position_sorted.bdg ~/Mus_musculus/mm10_2.genome /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position.bw

python /nethome/axy148/MACS/bin/macs2 callpeak -t /scratch/projects/bbc/BAM2BW/"$sample_name"_bam_sorted_by_position.bam -f BAM -g mm -n "$sample_name"_bam_mm --outdir /scratch/projects/bbc/Peak_chip_seq -B -q 0.01

EOF

#bsub -P Bioinformatics4count < ~/Script_bash/Run_"$sample_name"_to_bw.sh

bsub -e %J.err -P Bioinformatics4count < ~/Script_bash/Run_"$sample_name"_to_chip_seq_peak.sh

done < $1

#samtools sort 201348189-01.bam 201348189-01_bam_sorted_by_position
#genomeCoverageBed -ibam 201348189-01_bam_sorted_by_position.bam -bg -g mm10.genome > 201348189-01_bam_sorted_by_position.bdg
#bedtools slop -i 201348189-01_bam_sorted_by_position.bdg -g mm10.genome -b 0 
#bedGraphToBigWig 201348189-01_bam_sorted_by_position.bdg mm10_2.genome 201348189-01_bam_sorted_by_position_2.bw
#| bedClip stdin mm10.genome 201348189-01_bam_sorted_by_position.bdg.clip