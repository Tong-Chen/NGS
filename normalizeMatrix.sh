#!/bin/bash

#set -x

usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to compute Z-score of a matrix.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, the first column is the
 		colname, tab seperated)${bldred}[NECESSARY]${txtrst}
	-c	The number of columns you want to use
		[${txtred}Default 0 means all columns.${txtrst}]
	-r	The number of rows you want to use
		[${txtred}Default 0 means all rows.${txtrst}]
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
	-e	Execute or not[${bldred}Default TRUE${txtrst}]
EOF
}

file=
row=0
col=0
header='TRUE'
execute='TRUE'

while getopts "hf:r:c:z:e:" OPTION
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
			row=$OPTARG
			;;
		c)
			col=$OPTARG
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

cat <<END >${file}.r

data <- read.table(file="${file}", sep="\t", header=$header,
row.names=1, check.names=F)

data.m <- as.matrix(data)

if (${row} + ${col} == 0) {
	data.v <- as.vector(data.m)
	data.o <- (data-mean(data.v))/sd(data.v)
}else if (${row} == 0){
	col_name_bak <- colnames(data.m)
	data.c <- data.m[, 1:${col}]
	data.v <- as.vector(data.c)
	data.vo <- (data.c-mean(data.v))/sd(data.v)
	data.o <- cbind(data.vo, data[,-(1:${col})])
	colnames(data.o) <- col_name_bak
}else if (${col} == 0){
	row_name_bak <- rownames(data.m)
	data.r <- data.m[1:${row},]
	data.v <- as.vector(data.r)
	data.vo <- (data.r-mean(data.v))/sd(data.v)
	print(mean(data.v))
	print(sd(data.v))
	data.o <- rbind(data.vo, data[-(1:${row}),])
	rownames(data.o) <- row_name_bak
}

write.table(data.o, file="${file}.Zscore", 
row.names=T, col.names=NA, sep="\t", quote=F)


END

if [ "$execute" == "TRUE" ]; then
	Rscript ${file}.r
	sed '1s/^\t/label\t/' -i ${file}.Zscore
	rm -f ${file}.r
fi
