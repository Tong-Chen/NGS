#!/bin/bash

set -x
set -e
set -u

usage()
{
cat <<EOF >&2
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to download sra files from a file containing SRA
accession numbers and transfer SRA to fastq format.

The format of input file:

SRR2952884      Y25_root
SRR2952883      Y18_root
SRR2952882      Y12_root
SRR2952881      Y05_root

The second column will be treated as the prefix of final fastq files.

${txtbld}OPTIONS${txtrst}:
	-f	Data file with format described above${bldred}[NECESSARY]${txtrst}
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
EOF
}

file=

while getopts "hf:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

if [ -z $file ]; then
	usage
	exit 1
fi


#IFS="\t"

cat $file | while read -r -a array 
do
	sra="${array[0]}"
	name="${array[1]}"
	#echo $sra, $name
	#prefetch -v $sra
	#prefetch -v $sra
	#prefetch -v $sra
	#/bin/cp ~/ncbi/public/sra/${sra}.sra ${name}.sra
	while true
	do
		fastq-dump -v --split-3 --gzip ${sra}
		a=$?
		if [ "$a" == "0" ]; then break; fi
		sleep 5m
	done
	
	#/bin/cp ~/ncbi/public/sra/${sra}* .
	if [ "$a" == "0" ]
	then
		rename "${sra}" "${name}" ${sra}*
		/bin/rm ~/ncbi/public/sra/${sra}.sra
	fi
done
