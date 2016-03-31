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

This script is used to do word cloud plot.

The format of file given to -f:

>>Term	Freq	Others
>>first sentence	3	Others
>>word	4	Others
>>third	5	Others

For this file, -w should be <Term> and -F should be <Freq>. Other
columns will be ignored.



${txtbld}OPTIONS${txtrst}:
	-f	Data file ${bldred}[NECESSARY]${txtrst}
	-w	The name of word column.${bldred}[NECESSARY]${txtrst}
	-F	The name of frequency column.${bldred}[NECESSARY]${txtrst}
	-m	The minimum freq of words to be output.${bldred}[Default 1.3]${txtrst}
	-M	The maximum numer of words to be output.${bldred}[Default 120]${txtrst}
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
	-e	Execute the secript or not[${bldred}Default TRUE${txtrst}]

Eg: $0 -f at_leat_two_colum_file_with_header_line -w Term -F freq

EOF


}

file=
header='TRUE'
word=
freq=
execute='TRUE'
min_freq=1.3
max_count=120

while getopts "he:f:m:M:w:F:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		F)
			freq=$OPTARG
			;;
		m)
			min_freq=$OPTARG
			;;
		M)
			max_count=$OPTARG
			;;
		w)
			word=$OPTARG
			;;
		e)
			execute=$OPTARG
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

mid=".wordcloud"

cat <<EOF >${file}${mid}".r"

library(wordcloud)
library(RColorBrewer)

pal <- brewer.pal(9, "BuGn")
pal <- pal[-(1:6)]

data <- read.table(file="${file}", header=${header}, sep="\t",
quote="", comment="")

data\$${freq}[data\$${freq}==Inf] <- max(data\$${freq}[data\$${freq}!=Inf])

png("${file}${mid}.png", pointsize=9)
wordcloud(data\$${word}, data\$${freq}, scale=c(2,.9),
	min.freq=${min_freq}, max.words=${max_count}, random.order=T, 
	random.color=F, rot.per=0,
	colors=pal, ordered.colors=F,  use.r.layout=F,  fixed.asp=T)
dev.off()

EOF

if test "${execute}" == "TRUE"; then
	Rscript ${file}${mid}.r
fi
