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

echo "Begin ${op}---`date`" >>${op}.log
echo "Begin mapping---`date`" >>${op}.log
bwa aln -t 6 ${index} ${fq} >${op}.sai 2>>${op}.log
bwa samse ${index} ${op}.sai ${fq} > ${op}.sam 2>>${op}.log
/bin/rm -f ${op}.sai &
echo "End mapping---`date`" >>${op}.log
echo "Begin filter unmapped reads and sort---`date`" >>${op}.log
head -n 1000 ${op}.sam | grep '^@' >${op}.uniq.sam
grep 'XT:A:U' ${op}.sam >>${op}.uniq.sam
samtools view -F 4 -ubS ${op}.uniq.sam \
	| samtools sort - ${op}.uniq.map.sort 2>>${op}.log
/bin/rm -f ${op}.uniq.sam &
echo "Total reads:`samtools view -cS ${op}.sam`" >>${op}.sta
/bin/rm -f ${op}.sam &
echo "Total unique redundant reads:`samtools view -c ${op}.uniq.map.sort.bam`" >>${op}.sta
touch ${op}
pythonmail2.py ${op} ${op}
