#!/bin/bash

#set -x

usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to generate correlation coefficient matrix for 
all rows in one file.  

The input file is like:

	gene_label	samp1	sampl2	sampl3
	1	0.70845836	0.38921956	-1.4000962
	2	-1.54789160	0.92229273	-0.5476238
	3	0.07930535	-0.69559018	0.7442455
	4	0.50076958	1.23827656	0.8005727


If -t is FALSE, the output is:

	label	1           2          3          4
	1  1.00000000 -0.03157877 -0.7592889 -0.0333722
	2 -0.03157877  1.00000000 -0.6264518  0.9999984
	3 -0.75928887 -0.62645177  1.0000000 -0.6250521
	4 -0.03337220  0.99999839 -0.6250521  1.0000000

If -t is TRUE, the output is:

Var1	Var2	value
1	1	1.00000000
2	1	-0.03157877
3	1	-0.75928887
4	1	-0.03337220
1	2	-0.03157877
2	2	1.00000000
3	2	-0.62645177
4	2	0.99999839
1	3	-0.75928887
2	3	-0.62645177
3	3	1.00000000
4	3	-0.62505211
1	4	-0.03337220
2	4	0.99999839
3	4	-0.62505211
4	4	1.00000000


The parameters for logical variable are either TRUE or FALSE.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, the first column is the
 		rowname, all columns are tab seperated)${bldred}[NECESSARY]${txtrst}
	-m	The method.[${bldred}Default pearson, accept "kendall",
		"spearman".${txtrst}]
	-t	Transfer correlation coefficient matrix to three column file
		in format like below
		gene1	gene2	correlation_coefficient
		A	A	1
		A	B	0.9
		A	C	0.5
		A	D	0.2
	-e	Execute or not[${bldred}Default TRUE${txtrst}]
	-i	Install the required packages[${bldred}Default FALSE${txtrst}]
EOF
}

file=''
execute='TRUE'
ist='FALSE'
method='pearson'
transfer='FALSE'

while getopts "hf:t:m:e:i:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		t)
			transfer=$OPTARG
			;;
		m)
			method=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		i)
			ist=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done

midname=".cor"

if [ -z $file ] ; then
	echo 1>&2 "Please give filename."
	usage
	exit 1
fi


cat <<END >${file}${midname}.r

data <- read.table(file="$file", sep="\t", header=T, row.names=1)
attach(data)

cc <- t(cor(t(data), use="everything", method="$method"))

if ($transfer){

	if(${ist}){
		install.packages("reshape2", repo="http://cran.us.r-project.org")
	}
	library(reshape2)
	cc.m <- melt(cc)
	write.table(cc.m, file="${file}.$method.cor", sep="\t",
	col.names=T, row.names=F, quote=F)
}else{
	write.table(cc, file="${file}.$method.cor", sep="\t", col.names=NA,
	row.names=T, quote=F)
}

END

if [ "$execute" == "TRUE" ]; then
	Rscript ${file}${midname}.r
	sed -i 's/^\t/label\t/' ${file}.$method.cor
fi

#convert -density 200 -flatten ${file}${midname}.eps ${first}${midname}.png
