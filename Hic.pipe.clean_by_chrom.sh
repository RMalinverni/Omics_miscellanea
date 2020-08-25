#!/bin/bash


COR=8
JTDIR="~/soft/juicer_tools_1.22.01.jar" # defalut (in my computer)
GENOME=hg19

while getopts i:t:o:c:n:j: option
do
	case "${option}"
	in

	i) BAM=${OPTARG};;
	t) COR=${OPTARG};;
	o) OUTDIR=${OPTARG};;
	c) CHR=${OPTARG};;
	n) NAME=${OPTARG};;
  j) JTDIR=${OPTARG};;
  g) GENOME=${OPTARG};;

esac

done

if [ ! -d $OUTDIR ]; then
	
	mkdir $OUTDIR
	mkdir $OUTDIR/samtmp

else
	echo $OUTDIR already exist	

fi	

if [ ! -f $OUTDIR/$NAME.sort.bam ]; then
echo 'starting sorting for position of file:' $BAM
#samtools sort -@ $COR $BAM -o $NAME.sort.bam 
sambamba sort -t $COR -m 32G --tmpdir $OUTDIR/samtmp $BAM -o $NAME.sort.bam
echo 'starting indexing of file:' $BAM
#samtools index -@ $COR $NAME.sort.bam
sambamba index -n $COR $NAME.sort.bam 
else
	echo 'find $NAME.sort.bam in $OUTDIR'
fi	

mv $NAME.sort.* $OUTDIR  #move .bam and .bai files

cd $OUTDIR

echo 'create header:' $NAME.sort.bam

samtools view -H $NAME.sort.bam > header.txt
echo 'cleaning bam file for chr:' $CHR 
samtools view $NAME.sort.bam $CHR | awk '$7 ~ /=/' |cat header.txt - | samtools view -Sb - > $NAME.$CHR.clean.bam

samtools sort -@ $COR -n $NAME.$CHR.clean.bam -o $NAME.$CHR.bam 
echo 'starting hicup2juicer' 
~/soft/hicup_v0.7.3/Conversion/hicup2juicer  $NAME.$CHR.bam
echo 'starting juicer tool Pre  and creating .hic file:' $NAME.$CHR.hic 
java -Xmx62g -jar $JTDIR pre --threads $COR $NAME.$CHR.bam.prejuicer $NAME.$CHR.hic $GENOME
