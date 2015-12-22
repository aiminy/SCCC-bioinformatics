#Usage:sh ConvertBed2Bw_12_22_2015_2.sh input_bed_file reference_genome_name

#example:
#sh ~/Code/ConvertBed2Bw_12_22_2015_2.sh /media/DATA/ChrisWilliams/ChIPseq/eto2_peaks.bed mm10
#you will find outfiles like the folowing:
#/media/DATA/ChrisWilliams/ChIPseq/eto2_peaks.unsorted.bedGraph
#/media/DATA/ChrisWilliams/ChIPseq/eto2_peaks.sorted.bedGraph
#/media/DATA/ChrisWilliams/ChIPseq/eto2_peaks.bw

echo $1
sample_name=`echo "$1" | awk -F"." '{print $1}'`
echo "$sample_name"

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' $1 > "$sample_name".unsorted.bedGraph  
sort -k1,1 -k2,2n  "$sample_name".unsorted.bedGraph > "$sample_name".sorted.bedGraph
mkdir ReferenceGenomeChromSize
fetchChromSizes.sh $2 > ReferenceGenomeChromSize/"$2".genome.chrom.sizes
bedGraphToBigWig "$sample_name".sorted.bedGraph ReferenceGenomeChromSize/"$2".genome.chrom.sizes "$sample_name".bw
