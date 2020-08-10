#!/bin/bash

# Softwares requested
#
# fastqc	- https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# samtools	- http://samtools.sourceforge.net/ 
# STAR		- https://github.com/alexdobin/STAR
# Subread	- http://subread.sourceforge.net/
# Multiqc	- https://multiqc.info/


PROJECT_DIR="X"
INDEX="X"
GFF="X"
NAME_EXP="RNAseq_analysis"
CORES=4
PAIRED="paired"


while getopts d:n:g:i:c:p: option
	
do
	case "${option}"
in
	
	d) PROJECT_DIR=${OPTARG};;
	n) NAME_EXP=${OPTARG};;
	g) GFF=${OPTARG};;
	i) INDEX=${OPTARG};;
	c) CORES=${OPTARG};;
	p) PAIRED=${OPTARG};;


esac
done


if [ $GFF = "X" ]; then
	echo "you need to provide a valid GFF file path"
	exit
fi

if [ $INDEX = "X" ]; then
	echo "you need to provide a valid STAR Index file path"
	exit
fi

if [ $PROJECT_DIR = "X" ]; then
	echo "you need to provide a valid project directory file path"
	exit
fi


#  Create a folder structure

FASTQ_DIR=$PROJECT_DIR"/fastq"
BAM_DIR=$PROJECT_DIR"/bam/"
LOGS_DIR=$PROJECT_DIR"/logs/"
FC_DIR=$PROJECT_DIR"/FCounts/"
TABS_DIR=$PROJECT_DIR"/tabs/"
HTML_DIR=$PROJECT_DIR"/html/"
#GFF=/scratch/ijcres/References/Annotations/Hsapiens/GRCh38_annotation_gencode_v32/gencode.v32.annotation.gff3
NTRHEADS=$CORES
#INDEX=/scratch/ijcres/Indexes/STAR/GRCh38

mkdir $BAM_DIR
mkdir $FASTQ_DIR
mkdir $LOGS_DIR
mkdir $FC_DIR
mkdir $TABS_DIR
mkdir $HTML_DIR

echo "#####################################"
echo "name RNAseq experimet: " $NAME_EXP
echo "project directory: " $PROJECT_DIR
echo "type of experimet: " $PAIRED
echo "number of cores: " $CORES
echo "#####################################"

cd $PROJECT_DIR

#name of fastq need to be NAME_***.fastq.gz if single
#name of fastq need to be NAME_***1(2).fastq.gz if paired

if [ $PAIRED = "paired" ]; then
               		
	
	ALL_SAMPLE=$(ls *1.fastq.gz | cut -d "_" -f 1)

	echo samples names : $ALL_SAMPLE
	echo "#####################################"
	for  i in $(ls *_R1.fastq.gz);
	
	do 
		echo "file: "$i

		NAME=$(echo $i | cut -d '_' -f 1)

		echo "name: "$NAME
		fastqc $NAME"_R1.fastq.gz"
		fastqc $NAME"_R2.fastq.gz"
		
		zcat $NAME"_R1.fastq.gz" > f1.fastq
		zcat $NAME"_R2.fastq.gz" > f2.fastq
		

		STAR --runThreadN $NTRHEADS --genomeDir $INDEX --readFilesIn f1.fastq f2.fastq  \
		--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999  \
		--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
		--alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix $NAME \
		--outStd SAM | samtools view -bS - > $NAME"_paired.bam"

		# --genomeLoadLoadAndKeep ### use if use different STAR process ( you will upload only once)
	
		rm f1.fastq
		rm f2.fastq
		mv $NAME"_R1.fastq.gz" $FASTQ_DIR
		mv $NAME"_R2.fastq.gz" $FASTQ_DIR  
	done
	
	LIST_BAM=$(ls *.bam)
	featureCounts -a $GFF -p -B -C -T $CORES -o $NAME_EXP"_FCounts.csv" $LIST_BAM
fi

if [ $PAIRED = "single" ]; then
	


	ALL_SAMPLE=$(ls fastq.gz | cut -d "_" -f 1)
	echo "#####################################"
	for  i in $(ls *fastq.gz);
	
	do 
		echo "file: "$i

		NAME=$(echo $i | cut -d '_' -f 1)

		echo "name: "$NAME
		
		fastqc $i

		zcat $i > sample.fq


		STAR --runThreadN $NTRHEADS --genomeDir $INDEX --readFilesIn sample.fq  \
		--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999  \
		--outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
		--alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix $NAME \
		--outStd SAM | samtools view -bS - > $NAME"_single.bam"

		# --genomeLoadLoadAndKeep ### use if use different STAR process ( you 			will upload only once)
	
		rm sample.fq

		mv $NAME $FASTQ_DIR
		
	done
	
	LIST_BAM=$(ls *.bam)
	featureCounts -a $GFF  -B -C -T $CORES -o $NAME_EXP"_FCounts.csv" $LIST_BAM

fi



multiqc .

mv *.bam $BAM_DIR
mv *.csv $FC_DIR
mv *.summary $FC_DIR
mv *.html $HTML_DIR
mv *.tab $TABS_DIR
mv *.fastq.gz $FASTQ_DIR
mv *.*  $LOGS_DIR



