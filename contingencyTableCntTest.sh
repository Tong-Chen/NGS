#!/bin/bash

usage()
{
cat <<EOF
${txtcyn}
Usage:

fileformat(for statistical test from beginning, total 200 for each
column):
Gene	samp1	samp2
q	12	78
w	13	13
e	0	9
r	10	2

fileformat( for multiple test only, the name for p-value colum is <p>
by default. One may change to other variable by specifying the <-M>
parameter.)
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
	-c	Cols names. Must be accordant with each column and total
		amount of each column.
		[${bldred}Optional, format "'F161','F1MEF'" ${txtrst}]
	-s	Total number of each columns. 	
		[${bldred}NECESSARY, format "10000,10010" ${txtrst}]
	-m	Statistical test method. 
		[${bldred}Default Fisher test(fisher.test), other optional
		'hypergeometry (phyper)', 'chisq (chisq)', 'McNemar
		(incomplete)'${txtrst}]
	-n	Only execute multi-test. [Default FALSE, accept TRUE]
	-M	Specify the name of p-value column. [Default <p>, only used
		when -n is TRUE]	
	-a	The minimum allowed p-value for log process.[Default 1e-10
		which means the largest value for neg_log10_p is 10. This
		only affects the volcano plot (often gives a better
		volcano visularization)]
	-p	Multi-test adjustment method.
		[${bldred}Default BH, other optional 'BY',
		'bonferroni', 'hochberg', 'homel', 'no' means no multiple
		test.${txtrst}]
	-q	Add 1 for computing fold change.
		[${bldred}Default TRUE, accept a number like 0 (no adding) or other
		integer or float${txtrst}]
	-r	Remove rows contain zero only.
		[${bldred}Default TRUE,  recommeded for multi-test. Accept
		FALSE${txtrst}]
	-t	Parameters for select significant differences.
		[${bldred}Default "2,0.05,0.1" which means fold change no less than 2
		fold (Though 2 and 0.5 are equal,  only number larger than 1
		is accepted), p-value no larger than 0.05, fdr
		no larger than 0.1. Format "fc,p-calue,fdr")${txtrst}]
	-v	Volcano plot.[${bldred}Default FALSE, if TRUE invokes
		volcano.sh.${txtrst}]
	-P	Parameters for volcano plot except -x,-y,-s.[${bldred}Default
		empty${txtrst}]
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
staTest='fisher'
multi_only='FALSE'
multi='BH'
col_name=''
amount=''
add_number=1
removeZero='TRUE'
sig="2,0.05,0.1"
volcano='FALSE'
volcano_p=''
minimum="1e-10"
p_value='p'

while getopts "hf:k:c:s:a:m:M:n:p:q:r:t:v:P:z:e:" OPTION
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
		c)
			col_name=$OPTARG
			;;
		s)
			amount=$OPTARG
			;;
		a)
			minimum=$OPTARG
			;;
		m)
			staTest=$OPTARG
			;;
		M)
			p_value=$OPTARG
			;;
		n)
			multi_only=$OPTARG
			;;
		p)
			multi=$OPTARG
			;;
		q)
			add_number=$OPTARG
			;;
		r)
			removeZero=$OPTARG
			;;
		t)
			sig=$OPTARG
			;;
		v)
			volcano=$OPTARG
			;;
		P)
			volcano_p=$OPTARG
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

midname="."${staTest}


cat <<EOF >$file${midname}.r
data <- read.table("$file", header=$header,
sep="\t",row.names=1, comment.char="", check.names=${checkNames})

amount <- c($amount)
k <- amount[1] #first data column
t <- amount[2] #second data column

if (! ${multi_only}){

col_name <- c($col_name)

if ("$col_name" != "" ){
	colnames(data) <- col_name
}

#data <- as.matrix(data)

if ($removeZero == "TRUE") {
	data <- data[rowSums(data)!=0,]
}

if ("$staTest" == 'fisher' || "$staTest" == 'chisq'){
	p <- numeric(nrow(data))
	for(i in 1:length(p)) {
		q <- data[i,1]
		m <- data[i,2]
		p[i] <- ${staTest}.test(matrix(c(q,k-q,m,t-m), nrow=2,
		byrow=FALSE))\$p.value
		#if (i %% 1000 == 0){
		#	print(i)
		#}
	}
}else if ("$staTest" == 'phyper'){
	p <- phyper(data[,1], data[,2], t-data[,2],k)
}else {
	print("Unknown statistical method,  ${staTest}")
	quit()
} 

data\$p <- p
} #---end statistical test--------------------------
#---end statistical test--------------------------

time <- t / k

data\$fc <- (data[,1]+${add_number}) / (data[,2]+${add_number}) * time

if ("${multi}" != "no"){
	fdr <- p.adjust(data\$p, "${multi}")
	data\$fdr <- fdr
}



data\$log2_fc <- log(data\$fc) / log(2)
data\$neg_log10_p <- log(data\$p) / log(10) * (-1)

data\$neg_log10_p[data\$neg_log10_p<${minimum}] <- ${minimum}

sig <- c($sig)

if ("${multi}" != "no"){
	data\$sig <- (data\$fc <= 1/sig[1] | data\$fc >=sig[1]) & data\$p <=
		sig[2] & data\$fdr < sig[3]
} else{
	data\$sig <- (data\$fc <= 1/sig[1] | data\$fc >=sig[1]) & data\$p <=
		sig[2]
}
file="${file}${midname}.result"
write.table(data, file=file, sep="\t", col.names=NA, row.names=T, quote=F)
EOF

if [ "${execute}" = 'TRUE' ]; then
	Rscript $file${midname}.r
	if [ "${volcano}" = 'TRUE' ]; then
		volcano.sh -f ${file}${midname}.result -x log2_fc -y neg_log10_p \
		-s sig ${volcano_p}
	fi
fi
