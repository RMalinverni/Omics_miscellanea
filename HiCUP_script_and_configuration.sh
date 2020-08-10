#!/bin/bash

alias hicup='~/soft/hicup_v0.7.3/hicup'

OUTDIR=$(pwd)
CORES=8
QUIET=0
KEEP=0
ZIP=1
FORMAT=Sanger
LONG=800
SHORT=20
NAME="file"
RUN=1



while getopts F:f:o:c:n:r: option
        
do
        case "${option}"
in

        f) FQ1=${OPTARG};;
        F) FQ2=${OPTARG};;
        o) OUTDIR=${OPTARG};;
        c) CORES=${OPTARG};;
	n) NAME=${OPTARG};;
	r) RUN=${OPTARG};;
	bt) BT2=${OPTARG};;
	i) INDEX=${OPTARG};;
	d) DIGEST=${OPTARG};;
	fo) FORMAT=${OPTARG};;
	l) LONG=${OPTARG};;
	s) SHORT=${OPTARG};;
esac
done


if test [ ! -f $FQ1 ]; then
	echo fastq1 need to be provided
	exit
fi

if test [ ! -f $FQ2 ]; then
	echo fastq2 need to be provided
	exit
fi

if test [ ! -f $BT2 ]; then
	echo a path fo Bowtie2 need to be provided
	exit
fi

if test [ ! -f $INDEX ]; then
	echo a path fo Bowtie2 Index need to be provided
	exit
fi

if test [ ! -f $DIGEST ]; then
	echo a path fo Bowtie2 Index need to be provided
	exit
fi


mkdir $OUTDIR

cat > conf.$NAME.hicup.txt << EOF
#Configuration file for HiCUP $NAME

Outdir: $OUTDIR
Threads: $CORES
Quiet: $QUIET
Keep: $KEEP
Zip: $ZIP
Bowtie2: $BT2
Index: $INDEX
Digest: $DIGEST
Format: $FORMAT
Longest: $LONG
Shortest: $SHORT
######################################################
$FQ1
$FQ2
EOF


mv conf.$NAME.hicup.txt $OUTDIR

cd $OUTDIR
if [ $RUN = 1 ];
	then
		hicup --conf conf.$NAME.hicup.txt
	else
		echo 'selected the creation of conf file only, use (-r 1) to run the HICUP pipeline'
fi	       	




