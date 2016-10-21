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
	-f	Data file containing ENTREZ ids. Normally gene symbols would
		be OK too. ${bldred}[NECESSARY]${txtrst}
	-s	Species names as listed in 
		<http://www.genome.jp/kegg/catalog/org_list.html>.
		Like, hsa, oas, mmu, dme, ath.
		"anopheles","arabidopsis","bovine","canine","celegans",
		"chicken","chimp","coelicolor","ecolik12","ecsakai",
		"fly","gondii","human","malaria","mouse","pig","rat",
		"rhesus","xenopus","yeast" and "zebrafish".
	-p	P-value[${bldred}Default 0.05${txtrst}]
	-q	Q-value[${bldred}Default 0.2${txtrst}]
	-l	The minimum IDs required for performing KEGG analysis.
		[${bldred}Default 10.
		Enlarging this number may be needed when error information 
		<checkAtAssignment("logical", "ontology", "character"):
		‘ontology’ is not a slot in class "logical">
		${txtrst}]
	-r	Genome-wide annotation file like org.Hs.eg.db for hsa.
		Suitable to get gene symbols when entrez-IDs are used as
		input.  
		[${bldred}Optional, accept org.Hs.eg.db, org.At.tair.db,
		org.Ce.eg.db, org.Dr.eg.db, org.Dm.eg.db, org.Mm.eg.db${txtrst}]
	-i	Install required packages[${bldred}Default FALSE${txtrst}]
	-e	Execute the script[${bldred}Default TRUE${txtrst}]
EOF
}

file=
header='TRUE'
install='FALSE'
execute='TRUE'
p_value=0.05
q_value=0.2
species=
least_id=10
anno_db=CTCT

while getopts "hf:s:p:q:l:r:i:e:" OPTION
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
		r)
			anno_db=$OPTARG
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

cat <<END >${file}.clusterProfileKEGG.r

if (${install}) {
	#library(devtools)
	#devtools::install_github(c("GuangchuangYu/GOSemSim", "GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))
	source("http://bioconductor.org/biocLite.R")
	biocLite("BiocUpgrade") ## you may need this
	biocLite("clusterProfiler")
	if ("${anno_db}" != "CTCT"){
		biocLite("${anno_db}")
	}
}


if ("${anno_db}" != "CTCT"){
	library("${anno_db}")
} 

library(DOSE)
library(clusterProfiler)

	
data <- read.table("${file}", sep="\t", comment="", quote="")
colnames(data) <- c('gene', 'samp')
sampC <- unique(data\$samp)

for(samp in sampC) {
	id <- unique(data[data\$samp==samp, 1])
	id <- bitr_kegg(id, fromType='kegg', toType='uniprot', organism="${species}")[,2]
	if (length(id) < ${least_id}) {
		print(paste0("Less than 10 ids for ", samp, \
			". No KEGG enrichment performed."))
		output <- paste("${file}", samp, "KEGG", sep=".")
		system(paste0("touch ", output))
	} else {
		print(paste0("KEGG enrichment for ", samp))

		kk <- enrichKEGG(id, organism="${species}", keyType='uniprot', pvalueCutoff=${p_value},
			pAdjustMethod="BH", qvalueCutoff=${q_value})
		result <- setReadable(kk, "${anno_db}",  keytype="UNIPROT")
		#result <- summary(kk)
		output <- paste("${file}", samp, "KEGG.xls", sep=".")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=T)
	}
}


END

if [ "${execute}" == "TRUE" ]; then
	Rscript ${file}.clusterProfileKEGG.r
	afterRunclusterProfileKEGG.py -i ${file}
fi


