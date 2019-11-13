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

This script is used to do GO GSEA using clusterprofile.

The input file should be a two-column file with one header line, like
(system will sort ids by the second column in decresing order .
#------------------------
entrez_1	5    
entrez_2    4   
entrez_3    3  
entrez_4    -1
entrez_5    4
entrez_6    -2
entrez_7    2
entrez_8    1
#------------------------



${txtbld}OPTIONS${txtrst}:
	-f	Data file containing ENTREZ ids. ${bldred}[NECESSARY]${txtrst}
	-r	Genome-wide annotation file like org.Hs.eg.db for hsa.
		GO annotation was extracted from this annotation file.
		[${bldred}NECESSARY, accept org.Hs.eg.db (human),
		org.At.tair.db (arabidopsis), org.Rn.eg.db (rat), 
		org.Ce.eg.db (celegans), org.Dr.eg.db (zebrafish),
		org.Dm.eg.db (fly), org.Mm.eg.db (mouse)${txtrst}]
	-p	Adjusted P-value[${bldred}Default 0.1${txtrst}]
	-L	Loose FDR[${bldred}Default FALSE. ${txtrst}]
	-l	The minimum IDs required for performing GO GSEA.
		[${bldred}Default 0.
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

cat <<END >${file}.clusterProfileGogsea.r

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

gsea_plot_ct <- function(gseaResult, file, type){
	result_len <- length(gseaResult)
	for (i in 1:result_len) {
		id = gseaResult\$ID[i]
		description = gseaResult\$Description[i]
		tag = gsub(" ", "_", description)
		tmp = paste(file, type, tag, ".pdf", sep="_")
		pdf(tmp)
		gseaplot(gseaResult, id, title=description)
		dev.off()
	}
}

readable=TRUE

if ("${anno_db}" != "CTCT"){
	library("${anno_db}")
	readable=TRUE
} 

library(clusterProfiler)
	
data <- read.table("${file}", sep="\t", comment="", quote="")
data_g <- data[,2]
names(data_g ) <- data[,1]
data_g <- sort(data_g, decreasing=T)

typeL = c("BP", "MF", "CC")

if (length(data_g) < ${least_id}) {
	print(paste0("Less than ", ${least_id}, "ids for ", ". No GO GSEA performed."))
	for (type in typeL) {
		output <- paste("${file}", paste0(type, "_GO_gsea.xls"), sep=".")
		enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_DE_genes_or_less_than_5_de_genes;No_DE_genes;0/0;0/0;1;1;1;No_DE_genes;", length(data_g))
		result <- read.table(text=enrich, sep=";")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=F)
	}
} else {
	for (type in typeL) {
		if (${loose}) {
			MF <- gseGO(data_g, ont=type, "${anno_db}", pvalueCutoff=10,
				pAdjustMethod="BH")
			if(is.null(MF)){
				enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(data_g))
				result <- read.table(text=enrich, sep=";", header=T)
			} else {
				gsea_plot_ct(MF, "${file}", type)
				MF <- as.data.frame(MF) 
				MF <- MF[MF\$pvalue<=0.1, ]
				MF\$p.adjust <- p.adjust(MF\$pvalue,  method="BH")
				result <- MF[MF\$pvalue<${p_value} & MF\$p.adjust<${q_value}, ]
			}
		} else {
			MF <- gseGO(data_g, ont=type, "${anno_db}", pvalueCutoff=${p_value},
				pAdjustMethod="BH")
			if( dim(as.data.frame(MF))[1] ==0){
				enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(data_g))
				result <- read.table(text=enrich, sep=";", header=T)
			} else {
				gsea_plot_ct(MF, "${file}", type)
				result <- as.data.frame(MF)
			}
		}

		if( dim(as.data.frame(MF))[1] ==0){
			enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(data_g))
			result <- read.table(text=enrich, sep=";", header=T)
		}

		output <- paste("${file}", paste0(type, "_Gogsea.xls"), sep=".")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=T)
	}
}


END

if [ "${execute}" == "TRUE" ]; then
	Rscript ${file}.clusterProfileGogsea.r
	afterRunclusterProfileGogsea.py -i ${file}
fi


