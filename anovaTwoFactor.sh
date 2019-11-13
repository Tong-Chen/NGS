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

This script is used to do ANOVA analysis with two factors.

It accepts several types of input data.

1. Compute if there is significant difference for every compound in each two groups.

Data matrix

Compound	Ita1	Ita2	Chi1	Chi2	Ger1	Ger2 ...
pig	0	0	800	790.66089	200	200.49341 ...
rabbit	283.05845	145.22119	1224.42165	501.76551	1176.32261	497.01935 ...
monkey	5316.20028	5324.68285	8223.78233	6077.86028	10167.074	15420.62318 ...

Group file (two column file, header is needed but their real name will be omitted)

Individual	genotyep	treament
Ita1	Ita	0h
Ita2	Ita	0h
Ita3	Ita	0h
Ita4	Ita	1h
Ita5	Ita	1h
Ita6	Ita	1h
Chi1	Chi	0h
Chi2	Chi	0h
Chi3	Chi	0h
Chi4	Chi	1h
Chi5	Chi	1h
Chi6	Chi	1h
Ger1	Ger	0h
Ger2	Ger	0h
Ger3	Ger	0h
Ger4	Ger	1h
Ger5	Ger	1h
Ger6	Ger	1h


The output would be


Compound	genotype	treatment	genotype:treatment
monkey	0.177434068907696	0.733945335047997	0.0951758177204322
pig	0	0	0
rabbit	0.997701154807196	0.37525039314285	0.397819955522017


${txtbld}OPTIONS${txtrst}:
	-f	Data matrix ${bldred}[NECESSARY]${txtrst}
	-F	Filter data by mad value (Median Absolute Deviation) 
		${bldred}[Default 0 meaning no filter]${txtrst}
	-g	Group file. Normally the first column of group file
   		corresponds with the first column of data matrix to
		specify the group information (order not matter).	
		${bldred}[NECESSARY]${txtrst}
	-a	Factor 1, one coulmn name of group file ${bldred}[NECESSARY]${txtrst}
	-b	Factor 2, one column name of group file ${bldred}[NECESSARY]${txtrst}
	-Q	FDR filter. ${bldred}[Default 0.1]${txtrst}
	-I	(uppercase I)Do these two factors have interactions.
		${bldred}[Default TRUE]${txtrst}
EOF
}

file=
group=
mad_threshold=0
header='TRUE'
factor1=''
factor2=''
interaction="TRUE"
fdr=0.1

while getopts "hf:g:a:b:I:F:Q:" OPTION
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
			mad_threshold=$OPTARG
			;;
		g)
			group=$OPTARG
			;;
		a)
			factor1=$OPTARG
			;;
		b)
			factor2=$OPTARG
			;;
		I)
			interaction=$OPTARG
			;;
		Q)
			fdr=$OPTARG
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

if (${mad_threshold} > 0){
	mad_cmp = apply(test, 1, mad)
	test <- test[mad_cmp>${mad_threshold}, ]
}

test\$id <- rownames(test)
test_m <- melt(test, id.vars="id")
grp_file <- read.table("${group}", header=T, sep="\t", row.names = 1, quote="")
#colnames(grp_file) <- c("group")
# Here we use variable as column name because this column name exists in test_m
grp_file\$variable <- rownames(grp_file)
merged_data <- merge(test_m, grp_file, by="variable", all.x=TRUE)

if (${interaction}) {
	all_result <- by(merged_data, INDICES=merged_data\$id, FUN=function(m) {
	  m <- droplevels(m)
	  ct_aov <- aov(value ~ ${factor1}+${factor2}+${factor1}:${factor2}, data=m)
	  #ct_hsd <- TukeyHSD(ct_aov)
	  #result <- as.data.frame(ct_hsd\$group[,4, drop=F])
	  #colnames(result) <- paste0(unique(m\$id), "_", colnames(result))
	  result <- as.data.frame(summary(ct_aov)[[1]][["Pr(>F)"]][1:3], row.names=c("${factor1}","${factor2}","${factor1}_${factor2}"))
	  colnames(result) <- unique(m\$id)
	  t(result)
	})
	result <- do.call(rbind, unname(all_result))
} else {
	all_result <- by(merged_data, INDICES=merged_data\$id, FUN=function(m) {
	  m <- droplevels(m)
	  ct_aov <- aov(value ~ ${factor1}+${factor2}, data=m)
	  #ct_hsd <- TukeyHSD(ct_aov)
	  #result <- as.data.frame(ct_hsd\$group[,4, drop=F])
	  #colnames(result) <- paste0(unique(m\$id), "_", colnames(result))
	  result <- as.data.frame(summary(ct_aov)[[1]][["Pr(>F)"]][1:3], row.names=c("${factor1}","${factor2}"))
	  colnames(result) <- unique(m\$id)
	  t(result)
	})
	result <- do.call(rbind, unname(all_result))

}

result_adjusted = apply(result, 2, p.adjust,"BH")
colnames(result_adjusted) <- paste(colnames(result_adjusted), 'FDR', sep="_")

result <- cbind(result, result_adjusted)
result <- as.data.frame(result)

result <- cbind(id=rownames(result), result)
result_m <- melt(result, id.vars=c('id'))

result_m_0_1 <- result_m[result_m\$value<${fdr}, ]
write.table(result_m_0_1, file="${file}.${factor1}.${factor2}.sig.xls", sep="\t", quote=F, row.names=F, col.names=F)

system("s-plot vennDiagram -f ${file}.${factor1}.${factor2}.sig.xls -a ${factor1}_FDR -b ${factor2}_FDR -c ${factor1}_${factor2}_FDR")

#result
write.table(result, file="${file}.${factor1}.${factor2}.xls", sep="\t", quote=F, row.names=F)
END

Rscript ${file}.annnovar.r
