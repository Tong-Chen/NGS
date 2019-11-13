#!/bin/bash

#set -x
set -e
set -u

usage()
{
cat <<EOF >&2
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to transfer data matrix to markdown table format.

${txtbld}OPTIONS${txtrst}:
	-f	Data file ${bldred}[NECESSARY]${txtrst}
	-r	First <n> rows ${bldred}[Default all rows]${txtrst}
	-c	First <n> columns ${bldred}[Default all columns]${txtrst}
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
	-n	First column as row names[${bldred}Default <NULL> to use numbers as row names. Accept <1> to treat first column as row names.${txtrst}]
EOF
}

file=
nrow=-1
ncol=0
header='TRUE'
row_names='NULL'

while getopts "hf:r:c:n:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		r)
			nrow=$OPTARG
			;;
		c)
			ncol=$OPTARG
			;;
		n)
			row_names=$OPTARG
			;;
		z)
			header=$OPTARG
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

cat <<END >${file}.markdown.r
data <- read.table("${file}", sep="\t", check.names=F, quote="", comment="", header=${header}, row.names=${row_names}, nrows=${nrow})

if (${ncol} > 0){
	data <- data[, 1:${ncol}, drop=F]
}

a = knitr::kable(data, format="markdown")

write(a, file="${file}.markdown.table")

END

Rscript ${file}.markdown.r
/bin/rm -f ${file}.markdown.r
