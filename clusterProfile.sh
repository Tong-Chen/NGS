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

This script is used to do KEGG enrichment using clusterprofile.

The input file should be a two column file with no header line, like
#------------------------
entrez_1       samp_a
entrez_2       samp_a
entrez_3       samp_a
entrez_4       samp_a
entrez_5       samp_b
entrez_6       samp_b
entrez_7       samp_c
entrez_8       samp_c
#------------------------



${txtbld}OPTIONS${txtrst}:
	-f	Data file containing ENTREZ ids ${bldred}[NECESSARY]${txtrst}
	-s	Species names as listed in 
		<http://www.genome.jp/kegg/catalog/org_list.html>.
		Like, hsa, oas, mmu, dme, ath.
	-p	P-value[${bldred}Default 0.05${txtrst}]
	-q	Q-value[${bldred}Default 0.1${txtrst}]
	-l	The minimum IDs required for performing KEGG analysis.
		[${bldred}Default 10.
		Enlarging this number may be needed when error information 
		<checkAtAssignment("logical", "ontology", "character"):
		‘ontology’ is not a slot in class "logical">
		${txtrst}]
	-i	Install required packages[${bldred}Default FALSE${txtrst}]
	-e	Execute the script[${bldred}Default TRUE${txtrst}]
EOF
}

file=
header='TRUE'
install='FALSE'
execute='TRUE'
p_value=0.05
q_value=0.1
species=
least_id=10

while getopts "hf:s:p:q:l:i:e:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		s)
			species=$OPTARG
			;;
		i)
			install=$OPTARG
			;;
		p)
			p_value=$OPTARG
			;;
		q)
			q_value=$OPTARG
			;;
		l)
			least_id=$OPTARG
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

cat <<END >${file}.clusterProfile.r

if (${install}) {
	devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))
}

library(DOSE)
library(clusterProfiler)

data <- read.table("${file}", sep="\t", comment="", quote="")
colnames(data) <- c('gene', 'samp')
sampC <- unique(data\$samp)

for(samp in sampC) {
	id <- data[data\$samp==samp, 1]
	if (length(id) < ${least_id}) {
		print(paste0("Less than 10 ids for ", samp, \
			". No KEGG enrichment performed."))
		output <- paste("${file}", samp, "KEGG", sep=".")
		system(paste0("touch ", output))
	} else {
		print(paste0("KEGG enrichment for ", samp))
		kk <- enrichKEGG(id, species="${species}", pvalueCutoff=${p_value},
			pAdjustMethod="BH", qvalueCutoff=${q_value})
		result <- summary(kk)
		output <- paste("${file}", samp, "KEGG", sep=".")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=T)
	}
}


END

if [ "${execute}" == "TRUE" ]; then
	Rscript ${file}.clusterProfile.r
fi
