#! /bin/bash -l

# Softwares requested 
#
# cutadapt	- https://cutadapt.readthedocs.io/en/stable/
# bowtie2	- http://bowtie-bio.sourceforge.net/bowtie2/index.shtmlcutadapt
# fastqc	- https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# samtools	- http://samtools.sourceforge.net/
# deeptools 	- https://deeptools.readthedocs.io/en/develop/
# multiqc 	- https://multiqc.info/


FASTQ1="none"
FASTQ2="none"
PAIRED="single"
INDEX="X"
JOB_DIR="·"
BOWTIE2_OPT="--very-sensitive  -k 2 -t"
	
while getopts d:c:a:i:b: option
	
do
	case "${option}"
	in

	d) JOB_DIR=${OPTARG};;
	c) CORES=${OPTARG};;
	a) ADAPTERS=${OPTARG};;
	i) INDEX=${OPTARG};;
	b) BOWTIE2_OPT=${OPTARG};;
	
	esac
done

					
	date

	
	echo "###############################"
	echo "Output folder: " $JOB_DIR
	echo "Number of cores used: " $CORES
	echo "adapter's file: " $ADAPTERS
	echo "Bowtie2 index: "$INDEX
	echo "###############################"	
	
	echo "making structure directories"
	JOB_DIR=$FOLDER
	BAM_DIR=$JOB_DIR/bam_files
	BIGWIG_DIR=$JOB_DIR/bigwig_files
	FASTQ_DIR=$JOB_DIR/fastq_files
	HTML_DIR=$JOB_DIR/fastQC_files
	LOGS_DIR=$JOB_DIR/logs
	
	mkdir $BAM_DIR
	mkdir $BIGWIG_DIR	
	mkdir $FASTQ_DIR
	mkdir $HTML_DIR

	cd $JOB_DIR
	
	########## FastQ cleaning ##########	
	
	for i in (ls *.fastq.gz)
		
	do		
		echo "start cleaning files $i"
		
		NAME=${i%*.fastq.gz}
		
		cutadapt --cores $CORES -m 24 -a "file:"$ADAPTERS $i > $NAME.cutadapt.fq.gz
		
		mv $i /fastq_files	
		mv $NAME.cutadapt.fq.gz $NAME.clean.fq.gz
		
		fastqc $NAME.clean.fq.gz
				
	done


	########## Alignment ############# 	
	
	if [ PAIRED="single" ]; then
		
		ALL_SAMPLE=$(ls fastq.gz | cut -d "_" -f 1)

		echo "########################################"
		echo $ALL_SAMPLES
		echo "########################################"
		
		for i in $(ls *.clean.fq.gz)
						
		do
				
			INPUT_BOWTIE=$i			
			NAME=${i%*.clean.fq.gz}
	
			bowtie2  -x $INDEX -U $INPUT_BOWTIE $BOWTIE2_OPT -p $CORES | samtools view -bS - > $NAME".bam"
	

		done
	fi
	
	if [ $PAIRED = "paired" ]; then

		ALL_SAMPLE=$(ls *1.fastq.gz | cut -d "_" -f 1)

		echo "########################################"		
		echo $ALL_SAMPLES
		echo "########################################"

		for i in $(ls *1.clean.fq.gz);
		
		do
			NAME=${i%*1.clean.fq.gz}		
			INPUT_BOWTIE1=$NAME"1.clean.fq.gz"
			INPUT_BOWTIE2=$NAME"2.clean.fq.gz"
			
			bowtie2  -x $INDEX -1 $INPUT_BOWTIE1 -2 $INPUT_BOWTIE2 $BOWTIE2_OPT -p $CORES | samtools view -bS - > $NAME".bam"
		done
	fi		
	
	
	############# Cleaning and sorting ####################	
		
	samtools view -H $NAME".bam" > header.sam
	samtools view -f4 $NAME".bam" | cat header.sam - | samtools view -b - > $NAME"_unmapped.bam"
	samtools view -F 4 $NAME".bam" | grep -v "XS:" | cat header.sam - | \
	samtools view -b - > unique.bam
	mv unique.bam $NAME".bam"	
	rm header.sam

	echo "indexing file: "$NAME.bam	
	
	samtools sort  -o $NAME.sort.bam -@ $CORES $NAME.bam
	rm $NAME.bam	
	samtools index $NAME.sort.bam
			
			
	echo "removing duplicated reads: "$NAME.sort.bam
			
	samtools sort -@ $CORES -n -o namesort.bam $NAME.sort.bam
	samtools fixmate -m namesort.bam fixmate.bam
	samtools sort -@ $CORES -o positionsort.bam fixmate.bam
	samtools markdup -r positionsort.bam $NAME.markdup.bam
	samtools index $NAME.markdup.bam

	rm namesort.bam
	rm fixmate.bam
	rm positionsort.bam

	fastqc $NAME.markdup.bam  # make optional
	
	echo "calculating coverage and creating bigWig file"
			
	bamCoverage --binSize 50 -b $NAME.markdup.bam -o $NAME.bw
			
			
	rm *.zip 	
	
	rm	$FOLDER"/fastq_files/"$NAME".adapt.fq"
	rm	$NAME".sort.bam"

	multiqc .	
	
	mv *.bam $BAM_DIR
	mv *.html $HTML_DIR
	mv *.fastq.gz $FASTQ_DIR
	mv *.*  $LOGS_DIR

	echo "################################## Finish ###########################################"


	

	
	
