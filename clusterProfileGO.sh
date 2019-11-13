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

This script is used to do GO enrichment using clusterprofile.

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
	-r	Genome-wide annotation file like org.Hs.eg.db for hsa.
		GO annotation was extracted from this annotation file.
		[${bldred}NECESSARY, accept org.Hs.eg.db (human),
		org.At.tair.db (arabidopsis), org.Rn.eg.db (rat), 
		org.Ce.eg.db (celegans), org.Dr.eg.db (zebrafish),
		org.Dm.eg.db (fly), org.Mm.eg.db (mouse)${txtrst}]
	-s	Species names as listed below. [Depleted]
		"anopheles","arabidopsis","bovine","canine","celegans",
		"chicken","chimp","coelicolor","ecolik12","ecsakai",
		"fly","gondii","human","malaria","mouse","pig","rat","rno",
		"rhesus","xenopus","yeast" and "zebrafish".
	-p	P-value[${bldred}Default 0.05${txtrst}]
	-q	Q-value[${bldred}Default 0.2${txtrst}]
	-L	Loose FDR[${bldred}Default FALSE. ${txtrst}]
	-l	The minimum IDs required for performing GO analysis.
		[${bldred}Default 5.
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
q_value=0.2
species=
least_id=5
anno_db=CTCT
loose='FALSE'

while getopts "hf:s:p:q:l:L:r:i:e:" OPTION
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
		L)
			loose=$OPTARG
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

if [ -z $file ] || [ -z $anno_db ]; then
	usage
	exit 1
fi

cat <<END >${file}.clusterProfileGO.r

if (${install}) {
	#library(devtools)
	#devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))
	source("http://bioconductor.org/biocLite.R")
	biocLite("BiocUpgrade") ## you may need this
	biocLite("clusterProfiler")
	if ("${anno_db}" != "CTCT"){
		biocLite("${anno_db}")
	}
}

readable=TRUE

if ("${anno_db}" != "CTCT"){
	library("${anno_db}")
	readable=TRUE
} 

library(clusterProfiler)

	
data <- read.table("${file}", sep="\t", comment="", quote="")
colnames(data) <- c('gene', 'samp')
sampC <- unique(data\$samp)
typeL = c("BP", "MF", "CC")

for(samp in sampC) {
	id <- unique(data[data\$samp==samp, 1])
	if (length(id) < ${least_id}) {
		print(paste0("Less than 5 ids for ", samp, \
			". No GO enrichment performed."))
		for (type in typeL) {
			output <- paste("${file}", samp, paste0(type, "_GO.xls"), sep=".")
			enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_DE_genes_or_less_than_5_de_genes;No_DE_genes;0/0;0/0;1;1;1;No_DE_genes;", length(id))
			result <- read.table(text=enrich, sep=";")
			write.table(result, file=output, quote=F, sep="\t", row.names=F,
			col.names=F)
		}
	} else {
		print(paste0("GO enrichment for ", samp))
		for (type in typeL) {
			if (${loose}) {
				MF <- enrichGO(id, "${anno_db}", pvalueCutoff=10,
					pAdjustMethod="BH", qvalueCutoff=10, ont=type, 
					readable=readable)
				if(is.null(MF)){
					enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(id))
					result <- read.table(text=enrich, sep=";", header=T)
				} else {
					MF <- MF[MF\$pvalue<=0.1, ]
					MF\$p.adjust <- p.adjust(MF\$pvalue,  method="BH")
					result <- MF[MF\$pvalue<${p_value} & MF\$p.adjust<${q_value}, ]
				}
			} else {
				MF <- enrichGO(id, "${anno_db}", pvalueCutoff=${p_value},
					pAdjustMethod="BH", qvalueCutoff=${q_value}, ont=type, 
					readable=readable)
				if(is.null(MF)){
					enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(id))
					result <- read.table(text=enrich, sep=";", header=T)
				} else {
					result <- simplify(MF, cutoff=0.7, by="p.adjust", select_fun=min)
				}
			}

			#pdf(paste("${file}", samp, "MF_GO.pdf", sep="."))
			#dotplot(MF, showCategory=30)
			#dev.off()
			#result <- summary(MF)
			enrichedTerm <- dim(result)[1]
			if (is.null(enrichedTerm) || enrichedTerm < 1){
				enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(id))
				result <- read.table(text=enrich, sep=";", header=T)
			}

			output <- paste("${file}", samp, paste0(type, "_GO.xls"), sep=".")
			write.table(result, file=output, quote=F, sep="\t", row.names=F,
			col.names=T)
		}

#		BP <- enrichGO(id, "${anno_db}", pvalueCutoff=${p_value},
#			pAdjustMethod="BH", qvalueCutoff=${q_value}, ont="BP", 
#			readable=readable)
#		result <- simplify(BP, cutoff=0.7, by="p.adjust", select_fun=min)
#		#pdf(paste("${file}", samp, "BP_GO.pdf", sep="."))
#		#dotplot(BP, showCategory=30)
#		#dev.off()
#		#result <- summary(BP)
#		output <- paste("${file}", samp, "BP_GO.xls", sep=".")
#		write.table(result, file=output, quote=F, sep="\t", row.names=F,
#		col.names=T)
#
#		CC <- enrichGO(id, "${anno_db}", pvalueCutoff=${p_value},
#			pAdjustMethod="BH", qvalueCutoff=${q_value}, ont="CC", 
#			readable=readable)
#		result <- simplify(CC, cutoff=0.7, by="p.adjust", select_fun=min)
#		#pdf(paste("${file}", samp, "CC_GO.pdf", sep="."))
#		#dotplot(CC, showCategory=30)
#		#dev.off()
#		#result <- summary(CC)
#		output <- paste("${file}", samp, "CC_GO.xls", sep=".")
#		write.table(result, file=output, quote=F, sep="\t", row.names=F,
#		col.names=T)
	}
}


END

if [ "${execute}" == "TRUE" ]; then
	Rscript ${file}.clusterProfileGO.r
	afterRunclusterProfileGO.py -i ${file}
fi


