#!/bin/bash

#set -x
set -e
#set -u

usage()
{
cat <<EOF >&2
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to do GO enrichment using clusterprofile.

The input ID file should be a two column file with no header line, like
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

The input GO file is also a two-column file

#-------------------------
ont	gene
KEGG_GLYCOLYSIS_GLUCONEOGENESIS	gene1
KEGG_GLYCOLYSIS_GLUCONEOGENESIS	gene2
KEGG_GLYCOLYSIS_GLUCONEOGENESIS	gene3
KEGG_GLYCOLYSIS	gene1
KEGG_GLYCOLYSIS	gene4 
KEGG_CYP	gene5
#-------------------------



${txtbld}OPTIONS${txtrst}:
	-f	Data file containing GENE IDs. Normally gene IDs should match
		ids in supplied GO file. ${bldred}[NECESSARY]${txtrst}
	-g	GeneOntology or other type annotation file with format specified above. 
		${bldred}[NECESSARY]${txtrst}
	-t	Tag of this analysis such as MF, BP, CC, KEGG
		${bldred}[NECESSARY, no blank allowed]${txtrst}
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
anno=
header='TRUE'
install='FALSE'
execute='TRUE'
p_value=0.05
q_value=0.2
species=
least_id=5
anno_db=CTCT
loose='FALSE'

while getopts "hf:s:g:t:p:q:l:L:r:i:e:" OPTION
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

cat <<END >${file}.clusterProfile${tag}.r

if (${install}) {
	source("http://bioconductor.org/biocLite.R")
	biocLite("BiocUpgrade") ## you may need this
	biocLite("clusterProfiler")
}


library(clusterProfiler)

	
data <- read.table("${file}", sep="\t", comment="", quote="")
colnames(data) <- c('gene', 'samp')
sampC <- unique(data\$samp)

anno <- read.table("${anno}", sep="\t", comment="", quote="", header=T)
colnames(anno) <- c("ont", "gene")

for(samp in sampC) {
	id <- unique(data[data\$samp==samp, 1])
	if (length(id) < ${least_id}) {
		print(paste0("Less than 5 ids for ", samp, \
			". No GO enrichment performed."))
		output <- paste("${file}", samp, "${tag}.xls", sep=".")
		enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_DE_genes_or_less_than_5_de_genes;No_DE_genes;0/0;0/0;1;1;1;No_DE_genes;", length(id))
		result <- read.table(text=enrich, sep=";")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=F)
	} else {
		print(paste0("GO enrichment for ", samp))
		if (${loose}) {
			MF <- enricher(id, TERM2GENE=anno, pvalueCutoff=10,
				minGSSize=${least_id}, maxGSSize=2000, 
				pAdjustMethod="BH", qvalueCutoff=10)
			if(is.null(MF)){
				enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(id))
				result <- read.table(text=enrich, sep=";", header=T)
			} else {
				#MF <- simplify(MF, cutoff=0.7, by="p.adjust", select_fun=min)
				MF <- MF[MF\$pvalue<=0.1, ]
				MF\$p.adjust <- p.adjust(MF\$pvalue,  method="BH")
				result <- MF[MF\$pvalue<${p_value} & MF\$p.adjust<${q_value}, ]
			}
		} else {
			MF <- enricher(id, TERM2GENE=anno, pvalueCutoff=${p_value},
				minGSSize=${least_id}, maxGSSize=2000, 
				pAdjustMethod="BH", qvalueCutoff=${q_value})
			if(is.null(MF)){
				enrich=paste0("ID;Description;GeneRatio;BgRatio;pvalue;p.adjust;qvalue;geneID;Count\\nNo_enrichment;No_enrichment;0/0;0/0;1;1;1;No_enrichment;", length(id))
				result <- read.table(text=enrich, sep=";", header=T)
			} else {
				result <- MF
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

		output <- paste("${file}", samp, "${tag}.xls", sep=".")
		write.table(result, file=output, quote=F, sep="\t", row.names=F,
		col.names=T)

	}
}


END

if [ "${execute}" == "TRUE" ]; then
	Rscript ${file}.clusterProfile${tag}.r
	afterRunclusterProfileSelf.py -i ${file} -t "${tag}"
fi


