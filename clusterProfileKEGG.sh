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
	-L	Losse q-value computation. ${bldred}[Default FALSE]${txtrst}
	-l	The minimum IDs required for performing KEGG analysis.
		[${bldred}Default 10.
		Enlarging this number may be needed when error information 
		<checkAtAssignment("logical", "ontology", "character"):
		‘ontology’ is not a slot in class "logical">
		${txtrst}]
	-r	Genome-wide annotation file like org.Hs.eg.db for hsa.
		Suitable to get gene symbols when entrez-IDs are used as
		input.  
		[${bldred}Optional, accept org.Hs.eg.db, org.At.tair.db, org.Rn.eg.db, 
		org.Ce.eg.db, org.Dr.eg.db, org.Dm.eg.db, org.Mm.eg.db${txtrst}]
	-i	Install required packages[${bldred}Default FALSE${txtrst}]
	-e	Execute the script[${bldred}Default TRUE${txtrst}]
EOF
}

file=
header='TRUE'
install='FALSE'
loose='FALSE'
execute='TRUE'
p_value=0.05
q_value=0.2
species=
least_id=10
anno_db=CTCT

while getopts "hf:s:p:L:q:l:r:i:e:" OPTION
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
		L)
			loose=$OPTARG
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
	# For lateset version
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
	if ("${species}" == "ath") {
		id <- bitr_kegg(id, fromType='ncbi-geneid', toType='kegg', organism="${species}")[,2]
	}# else {
	#	id <- bitr_kegg(id, fromType='kegg', toType='uniprot', organism="${species}")[,2]
	#}
	#id <- bitr_kegg(id, fromType='ncbi-geneid', toType='uniprot', organism="${species}")[,2]
	#id <- bitr_kegg(id, fromType='ncbi-geneid', toType='kegg', organism="${species}")[,2]
	if (length(id) < ${least_id}) {
		print(paste0("Less than 10 ids for ", samp, \
			". No KEGG enrichment performed."))
		output <- paste("${file}", samp, "KEGG.xls", sep=".")
		enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_DE_genes_or_less_than_5_de_genes;No_DE_genes;0/0;0/0;1;1;1;No_DE_genes;", length(id))
		result <- read.table(text=enrich, sep=";")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=F)
		#output <- paste("${file}", samp, "KEGG", sep=".")
		#system(paste0("touch ", output))
	} else {
		print(paste0("KEGG enrichment for ", samp))
		
		if (${loose}) {
			#kk <- enrichKEGG(id, organism="${species}", keyType='ncbi-geneid', pvalueCutoff=0.1,
			#kk <- enrichKEGG(id, organism="${species}", keyType='uniprot', pvalueCutoff=0.1,
			kk <- enrichKEGG(id, organism="${species}", keyType='kegg', pvalueCutoff=0.1,
				pAdjustMethod="BH", qvalueCutoff=10)
		} else {
			#kk <- enrichKEGG(id, organism="${species}", keyType='ncbi-geneid', pvalueCutoff=${p_value},
			#kk <- enrichKEGG(id, organism="${species}", keyType='uniprot', pvalueCutoff=${p_value},
			kk <- enrichKEGG(id, organism="${species}", keyType='kegg', pvalueCutoff=${p_value},
				pAdjustMethod="BH", qvalueCutoff=${q_value})
		}
		if (is.null(kk)) {
			enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(id))
			result <- read.table(text=enrich, sep=";", header=T)
		} else {
			if ("${species}" == "ath") {
				result <- setReadable(kk, "${anno_db}", keytype="TAIR")
				#result <- kk
			} else {
				result <- setReadable(kk, "${anno_db}", keytype="ENTREZID")
				#result <- kk
			}

			if (${loose}) {
				result <- result[result\$pvalue<=0.1, ]
				result\$p.adjust <- p.adjust(result\$pvalue,  method="BH")
				result <- result[result\$pvalue<${p_value} & result\$p.adjust<${q_value}, ]
			}
			#result <- summary(kk)

			enrichedTerm <- dim(result)[1]
			if (is.null(enrichedTerm) || enrichedTerm < 1){
				enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(id))
				result <- read.table(text=enrich, sep=";", header=T)
			}
		}
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


