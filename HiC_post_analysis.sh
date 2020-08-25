#!/bin/sh


DIRPIPE=$(pwd)
CHRS='chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'
CORES=8
MOD1='no'
MOD2='no'
MOD3='no'
MOD4='no'
MOD5='no'
JTDIR='~/soft/juicer_tools_1.22.01.jar'

while getopts i:t:o:c:n:m1:m2:m3:m4:m5:j:d: option
do
	case "${option}"
	in

	i) ENCBAM=${OPTARG};;
	t) COR=${OPTARG};;
	o) OUTDIR=${OPTARG};;
	c) CHRS=${OPTARG};;
	n) NAME=${OPTARG};;
  m1) MOD1=${OPTARG};;
  m2) MOD2=${OPTARG};;
  m3) MOD3=${OPTARG};;
  m4) MOD4=${OPTARG};;
  m5) MOD5=${OPTARG};;
  j) JTDIR=${OPTARG};;
  d) DIRPIPE=${OPTARG};;

esac

done


alias juicertools='java -Xmx32g -jar $JTDIR'


#MODULE1 cleaning and make hic file for each chromosomes

if [ $MOD1 = 'yes' ];
	then
		for CHR in $CHRS 
		do
		sh  $DIRPIPE'/Hic.pipe.clean_by_chrom.sh' -i $ENCBAM -t $CORES -o $OUTDIR -c $CHR -n $NAME
		done	
	else
		echo 'Module 1 (Cleaning bam and matrix construction) is skipped'
fi
		
cd $OUTDIR
	
mkdir hic
mkdir bam
mkdir prejuicer
	
mv *.hic hic
mv *.prejuicer prejuicer
mv *.bam bam

cd hic

#MODULE2 TAD calling using 2 different parameters
if [ $MOD2 = 'yes' ];
	then
	echo '########### Activate Module 2 TADs call ############'
	for CHR in $CHRS
		do			
			echo 'TAD calling on: ' $CHR
			juicertools arrowhead -c $CHR -m 2000 -r 10000 -k KR --threads $CORES $NAME.$CHR.hic  $NAME.$CHR.out
			juicertools arrowhead -c $CHR -m 2000 -r 5000 -k KR --threads $CORES $NAME.$CHR.hic  $NAME.$CHR.out		
		done
	else
		echo 'Module 2 (TAD calling) is skipped'

fi

#MODULE 3 Loop calling using 3 different parameters		
if [ $MOD3 = 'yes' ];
	then
	echo '########### Activate Module 3 loops call ############'	
	for CHR in $CHRS
		do
			echo 'Loop calling on: ' $CHR
			juicertools hiccups -c $CHR -m 1024 -k KR -r 5000,10000,25000 -f .1,.1,.1 -p 4,2,1 -i 7,5,3 -t 0.02,1.5,1.75,2 -d 20000,20000,40000 --threads $CORES $NAME.$CHR.hic $NAME.$CHR.hiccups_loops
		done
	else
		echo 'Module 3 (Loop calling) is skipped'
fi		

if [ $MOD4 = 'yes' ];
	then
	echo '########### Activate Module 4 compartments call ############' 	
	for CHR in $CHRS
		do
			echo 'Compartments call on:' $CHR
			juicertools eigenvector KR $NAME.$CHR.hic $CHR BP 1000000 $NAME.$CHR.compartments_1M.txt
			juicertools eigenvector KR $NAME.$CHR.hic $CHR BP 500000 $NAME.$CHR.compartments_500K.txt
			juicertools eigenvector KR $NAME.$CHR.hic $CHR BP 100000 $NAME.$CHR.compartments_100K.txt
		done
	else
		echo 'Module 4 (Compartments calling) is skipped'	

fi

if [ $MOD5 = 'yes' ];
	then
	echo '########### Activate Module 5 Dump hic files  ############'	
	for CHR in $CHRS
		do
			echo 'Dump file: '  $NAME.$CHR.hic
		       juicertools dump observed KR $NAME.$CHR.hic $CHR $CHR BP 10000  $NAME.$CHR.dump.txt
	       	done
	else
 		echo 'Module 5 (Dump files) is skipped'
fi			



mkdir ../TADs
mkdir ../loops
mkdir ../compartments
mkdir ../dumps

mv *hiccups* ../loops
mv *.out ../TADs
mv *compartments* ../compartments
mv *dump* ../dumps

echo '########### Done!  ############'
