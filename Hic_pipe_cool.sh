#!/bin/bash

# Software requested
#
# samtools	- http://samtools.sourceforge.net/
# pairtools 	- https://pairtools.readthedocs.io/en/latest/
# BWA 		- http://bio-bwa.sourceforge.net/
# cooler 	- https://github.com/mirnylab/cooler






JOB_DIR="."
CORES=4
NAME="HiC_test"

while getopts F:f:d:c:n:i:s: option
	
do
	case "${option}"
	in
	
	F) FQ1=${OPTARG};;
	f) FQ2=${OPTARG};;
	d) JOB_DIR=${OPTARG};;
	c) CORES=${OPTARG};;
	n) NAME=${OPTARG};;
	i) INDEX=${OPTARG};;
	s) CHR_SIZE=${OPTARG};;	

	esac
done


cd $JOB_DIR

if test [ ! -f $FQ1 ]; then
	echo fastq1 need to be present in $(pwd)
	exit
fi

if test [ ! -f $FQ2 ]; then
	echo fastq2 need to be present in $(pwd)
	exit
fi

if test [ ! -f $CHR_SIZE ]; then
	echo chromosome size need to be specified
	exit
fi


mkdir pairs
mkdir bam
mkdir cooler


echo aligment using BWA

bwa mem -SP5M -t $CORES $INDEX $FQ1 $FQ2 | samtools view -bhS - > $NAME".bam"

echo filtering 

cd $W_DIR

samtools view -h $NAME.bam |  pairtools parse -c $CHR_SIZE -o $NAME.pair.gz

pairtools sort --nproc 28 -o $NAME.sorted.pair.gz $NAME.pair.gz
pairtools dedup --mark-dups -o $NAME.dedup.pair.gz $NAME.sorted.pair.gz
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' -o $NAME.filtered.pair.gz $NAME.dedup.pair.gz
pairtools split --output-pairs $NAME.output.pairs.gz  $NAME.filtered.pair.gz


pairix -f $NAME.output.pairs.gz

rm $NAME.pair.gz
rm $NAME.sorted.pair.gz
rm $NAME.dedup.pair.gz 
rm $NAME.filtered.pair.gz


cooler cload pairix $CHR_SIZE:1000 $NAME.output.pairs.gz $NAME.cool
cooler balance $NAME.cool
cooler zoomify $NAME.cool

mv *.cool /cooler
mv *.pairs.gz /pairs
mv *.bam /bam


