#!/bin/bash

set -x
set -e
set -u

if test $# -lt 3; then
	echo 1>&2 "Usage $0 index output_prefix fastq"
	exit 1
fi

index=$1
op=$2
fq=$3

echo "Begin mapping ${op}---`date`" >>${op}.log
bowtie -n2 -l28 -e70 --chunkmbs 300 -a -m1 -S --best \
	--strata -p6 $(index) $(fq) ${op}.sam>${op}.log 2>&1 
echo "End mapping ${op}---`date`" >>${op}.log
samtools view -F 4 -ubS ${op}.sam | samtools sort - ${op}.uniq.map.sort
echo "Total reads:`samtools view -cS ${op}.sam`" >>${op}.sta
echo "Total unique redundant reads:`samtools view -c ${op}.uniq.map.sort.bam`" >>${op}.sta
/bin/rm -f ${op}.sam &
pythonmail2.py ${op} ${op}
