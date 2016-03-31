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
	-M	Specify the name of p-value column. [Default <FDR>]	
	-L	Specify the name of log p-value column. [Default <neg_log10FDR]	
	-p	Multi-test adjustment method.
		[${bldred}Default BH, other optional 'BY',
		'bonferroni', 'hochberg', 'homel', 'no' means no multiple
		test.${txtrst}]
	-o	Output file.
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
p_value='FDR'
logp="neg_log10FDR"
output=
while getopts "hf:k:c:s:a:L:m:M:n:o:p:q:r:t:v:P:z:e:" OPTION
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
		L)
			logp=$OPTARG
			;;
		n)
			multi_only=$OPTARG
			;;
		o)
			output=$OPTARG
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

midname=".fdr"


cat <<EOF >$file${midname}.r
data <- read.table("$file", header=$header,
sep="\t", quote="", comment.char="", check.names=${checkNames})

fdr <- p.adjust(data\$${p_value}, "${multi}")
data\$${p_value} <- fdr
data\$${logp} <- (-1) * log10(data\$${p_value})

file="${output}"
write.table(data, file=file, sep="\t", col.names=T, row.names=F, quote=F)
EOF

if [ "${execute}" = 'TRUE' ]; then
	Rscript $file${midname}.r
	/bin/rm -f $file${midname}.r
fi
