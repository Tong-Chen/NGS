#!/bin/bash

usage()
{
cat <<EOF
${txtcyn}
Usage:

fileformat( for multiple test only, the name for p-value colum must be 'p')
Gene	samp1	samp2	p
q	12	78	7.52661553061e-11
w	13	13	1.0
e	0	9	0.00356304943458
r	10	2	0.0357799098682

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to compute statistical parameter line by line by
supplied methods.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, the first column is the
 		rowname, tab seperated,no dup)${bldred}[NECESSARY]${txtrst}
	-k	If the names of your rows and columns startwith numeric value,
		this can be set to FALSE to avoid modifying these names to be
		illegal variable names. But duplicates can not be picked out.
		[${bldred}Default TRUE${txtrst}]
		Accept FALSE.
	-p	Multi-test adjustment method.
		[${bldred}Default BH, other optional 'BY',
		'bonferroni', 'hochberg', 'homel', 'no' means no multiple
		test.${txtrst}]
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
		Accept FALSE.
	-e	Execute or not[${bldred}Default TRUE${txtrst}]
		Accept FALSE.
EOF
}

file=
checkNames='TRUE'
header='TRUE'
execute='TRUE'
multi='BH'

while getopts "hf:k:p:z:e:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		k)
			checkNames=$OPTARG
			;;
		p)
			multi=$OPTARG
			;;
		z)
			header=$OPTARG
			;;
		e)
			execute=$OPTARG
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


cat <<EOF >$file.r
data <- read.table("$file", header=$header,
sep="\t", quote="", comment.char="", check.names=${checkNames})

fdr <- p.adjust(data\$p, "${multi}")
data\$fdr <- fdr

file="${file}.${multi}"
write.table(data, file=file, sep="\t", col.names=T, row.names=F, quote=F)
EOF

if [ "${execute}" = 'TRUE' ]; then
	Rscript $file${midname}.r
fi
