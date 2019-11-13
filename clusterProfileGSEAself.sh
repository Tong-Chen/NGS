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

The input file should be a two-column file without header lines, like
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

Annotation file should be a two column file

ont	gene
KEGG_GLYCOLYSIS_GLUCONEOGENESIS	gene1
KEGG_GLYCOLYSIS_GLUCONEOGENESIS	gene2
KEGG_GLYCOLYSIS_GLUCONEOGENESIS	gene3
KEGG_GLYCOLYSIS	gene1
KEGG_GLYCOLYSIS	gene4 
KEGG_CYP	gene5


${txtbld}OPTIONS${txtrst}:
	-f	Data file containing ENTREZ ids. ${bldred}[NECESSARY]${txtrst}
	-g	GeneOntology or other type annotation file with format specified above. 
		${bldred}[NECESSARY]${txtrst}
	-t	Tag of this analysis such as MF, BP, CC, KEGG
		${bldred}[NECESSARY, no blank allowed]${txtrst}
	-p	Adjusted P-value[${bldred}Default 0.2${txtrst}]
	-l	The minimum IDs required for performing GSEA.
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
anno=
tag=
header='TRUE'
install='FALSE'
execute='TRUE'
p_value=0.05
q_value=0.2
species=
least_id=5
loose='FALSE'

while getopts "hf:g:t:s:p:q:l:r:i:e:" OPTION
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
			anno=$OPTARG
			;;
		t)
			tag=$OPTARG
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

if [ -z $file ] || [ -z $anno ]; then
	usage
	exit 1
fi

cat <<END >${file}.clusterProfileGSEA${tag}.r

if (${install}) {
	source("http://bioconductor.org/biocLite.R")
	biocLite("BiocUpgrade") ## you may need this
	biocLite("clusterProfiler")
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


library(clusterProfiler)
	
data <- read.table("${file}", sep="\t", comment="", quote="", header=F, row.names=NULL)
data_g <- data[,2]
names(data_g ) <- data[,1]
data_g <- sort(data_g, decreasing=T)

anno <- read.table("${anno}", sep="\t", comment="", quote="", header=T, row.names=NULL)
colnames(anno) <- c("ont", "gene")

typeL <- c("${tag}")

if (length(data_g) < ${least_id}) {
	print(paste0("Less than ", ${least_id}, "ids for ", ". No  GSEA performed."))
	for (type in typeL) {
		output <- paste("${file}", paste0(type, "_gsea.xls"), sep=".")
		enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_DE_genes_or_less_than_5_de_genes;No_DE_genes;0/0;0/0;1;1;1;No_DE_genes;", length(data_g))
		result <- read.table(text=enrich, sep=";")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=F)
	}
} else {
	for (type in typeL) {
		if (${loose}) {
			MF <- GSEA(data_g, TERM2GENE=anno, pvalueCutoff=${p_value},
				pAdjustMethod="BH", exponent = 1, nPerm = 1000,  
				minGSSize = 10, maxGSSize = 5000,  
				verbose = TRUE, by = "fgsea")
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
			MF <- GSEA(data_g, TERM2GENE=anno, pvalueCutoff=${p_value},
				pAdjustMethod="BH", exponent = 1, nPerm = 1000,  
				minGSSize = 10, maxGSSize = 5000,  
				verbose = TRUE, by = "fgsea")
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

		output <- paste("${file}", paste0(type, "_gsea.xls"), sep=".")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=T)
	}
}


END

if [ "${execute}" == "TRUE" ]; then
	Rscript ${file}.clusterProfileGSEA${tag}.r
	#afterRunclusterProfileGogsea.py -i ${file}
fi


