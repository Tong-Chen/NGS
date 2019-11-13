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

This script is used to do ANOVA analysis.

It accepts several types of input data.

1. Compute if there is significant difference for every compound in each two groups.

Data matrix

Compound	Ita1	Ita2	Chi1	Chi2	Ger1	Ger2
pig	0	0	800	790.66089	200	200.49341
rabbit	283.05845	145.22119	1224.42165	501.76551	1176.32261	497.01935
monkey	5316.20028	5324.68285	8223.78233	6077.86028	10167.074	15420.62318

Group file (two column file, header is needed but their real name will be omitted)

Individual	group
Ita1	Ita
Ita2	Ita
Ita3	Ita
Chi1	Chi
Chi2	Chi
Chi3	Chi
Ger1	Ger
Ger2	Ger
Ger3	Ger

The output would be


Compound	Ger-Chi	Ita-Chi	Ita-Ger
monkey	0.177434068907696	0.733945335047997	0.0951758177204322
pig	0	0	0
rabbit	0.997701154807196	0.37525039314285	0.397819955522017


${txtbld}OPTIONS${txtrst}:
	-f	Data matrix ${bldred}[NECESSARY]${txtrst}
	-g	Group file. Normally the first column of group file
   		corresponds with the first column of data matrix to
		specify the group information (order not matter).	
		${bldred}[NECESSARY]${txtrst}
EOF
}

file=
group=
header='TRUE'

while getopts "hf:g:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		g)
			group=$OPTARG
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

cat <<END >${file}.annnovar.r
library(reshape2)
library(data.table, quietly=T)
test <- read.table("${file}", header=T, row.names=1, quote="")
test\$id <- rownames(test)
test_m <- melt(test, id.vars="id")
grp_file <- read.table("${group}", header=T, row.names = 1, quote="")
colnames(grp_file) <- c("group")
grp_file\$variable <- rownames(grp_file)
merged_data <- merge(test_m, grp_file, by="variable", all.x=TRUE)
all_result <- by(merged_data, INDICES=merged_data\$id, FUN=function(m) {
  m <- droplevels(m)
  ct_aov <- aov(value ~ group, data=m)
  ct_hsd <- TukeyHSD(ct_aov)
  result <- as.data.frame(ct_hsd\$group[,4, drop=F])
  #colnames(result) <- paste0(unique(m\$id), "_", colnames(result))
  colnames(result) <- unique(m\$id)
  t(result)
})
result <- do.call(rbind, unname(all_result))
#result
write.table(result, file="${file}.anovar.xls", sep="\t", quote=F)
END

Rscript ${file}.annnovar.r
