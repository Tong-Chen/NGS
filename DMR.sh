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

This script is used to find differentially methylated regions and
differentially methylated loci using <methylKit> and <edmr>.

1. sampleFile

#--The header lines should not be changes---
#--sam: bismark mapped SAM files
#--name: one unqiue name to represent this file, each replicate should 
#		 have its own name.
#--group: specify  the groups of each file. Normally replicates belong to
#		  the same group.

#*****FILE CONTENT****when -s is FALSE******************
#sam	name	group
#T_6/T_6.sort.rmdup.sam	T_6	0
#C_9/C_9.sort.rmdup.sam	C_9	1
#*****FILE CONTENT****when -s is FALSE******************

#*****FILE CONTENT****when -s is TRUE******************
#sam	name	group
#T_6/T_6_CpG.txt	T_6	0
#C_9/C_9_CpG.txt	C_9	1
#*****FILE CONTENT****when -s is TRUE******************

${txtbld}OPTIONS${txtrst}:
	-f	A sample file with content described above. ${bldred}[NECESSARY]${txtrst}
	-s	Using methyl call results instead of calling methylation from bismark.
		${bldred}[Default FALSE]${txtrst}
	-o	Specify the output folder.${bldred}[NECESSARY]${txtrst}
	-p	Specify the output prefix.${bldred}[NECESSARY]${txtrst}
	-a	Specify the genome assembl version.
		${bldred}[NECESSARY, like hg19, mm9, oar3.1]${txtrst}
	-c	Minimum read coverage to call a methylation status for a base.
		${bldred}[Default 10]${txtrst}
	-t	Type of analysis to do. 
		${bldred}[Default "'CpG','CHG','CHH'"]${txtrst}
	-i 	Install required packages.${bldred}[Default FALSE, accept TRUE]${txtrst}
EOF
}

file=
called='FALSE'
output_dir=
prefix=
min_read_cov=10
header='TRUE'
install='FALSE'
assembl=
typeL="'CpG','CHG','CHH'"

while getopts "ha:c:f:i:o:p:s:t:" OPTION
do
case $OPTION in
	h)
		usage
		exit 1
		;;
	a)
		assembl=$OPTARG
		;;
	c)
		min_read_cov=$OPTARG
		;;
	f)
		file=$OPTARG
		;;
	i)
		install=$OPTARG
		;;
	o)
		output_dir=$OPTARG
		;;
	p)
		prefix=$OPTARG
		;;
	s)
		called=$OPTARG
		;;
	t)
		typeL=$OPTARG
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

mid=".DMR"

mkdir -p ${output_dir}

cat <<END >${output_dir}/${prefix}${mid}.r

if (${install}){
	# https://github.com/ShengLi/edmr
	# https://github.com/al2na/methylKit 
	install.packages( c("data.table",  "mixtools",  "devtools"))
	source("http://bioconductor.org/biocLite.R")
	biocLite(c("GenomicRanges", "IRanges"))
	# install from github
	library(devtools)
	install_github("ShengLi/edmr", build_vignettes=FALSE)
	# install the development version from github
	install_github("al2na/methylKit", build_vignettes=FALSE)
}

library(methylKit)
library(edmr)

sampleFile <- read.table("${file}", header=T, sep="\t",quote="")
bamL <- as.list(as.vector(sampleFile\$sam))
name <- as.list(as.vector(sampleFile\$name))
group <- as.vector(sampleFile\$group)

typeL <- c(${typeL})

dmr = function(type, pbj, name){
	if (type == "CpG") {
		meth <- unite(pbj, destrand=FALSE)
	} else {
		meth <- unite(pbj, destrand=FALSE)
	}
	corre_plot <- paste0("${output_dir}/${prefix}", ".", type, ".correlation.png" )
	png(corre_plot)
	getCorrelation(meth, plot=T)
	dev.off()
	length_name <- length(name[[1]])
	for(i in 1:length_name){
		label <- name[[1]][i]

		methylation_plot <- paste0("${output_dir}/${prefix}",".",label,'.',type,"_methylation.png")
		png(methylation_plot)
		getMethylationStats(pbj[[i]], plot=T, both.strands=F)
		dev.off()

		coverage_plot <- paste0("${output_dir}/${prefix}", ".", label, '.', type, "_coverage.png")
		png(coverage_plot)
		getCoverageStats(pbj[[i]], plot=T, both.strands=F)
		dev.off()
	}

	myDiff <- calculateDiffMeth(meth,  num.cores=30)
	myDiff_ct <- myDiff
	colname_myDiff <- colnames(myDiff)
	colname_myDiff[7] <- "meth_diff" 
	colname_myDiff -> colnames(myDiff_ct)
	
	first <- name[[1]][1]
	second <- name[[1]][2]
	prefix <- paste0("${output_dir}/", "${prefix}.", first, ".vs.", second)

	all <- cbind(meth, pvalue=myDiff_ct\$pvalue,
	qvalue=myDiff_ct\$qvalue, meth_diff=myDiff_ct\$meth_diff)
	colnames(all) <- c("chr", "start", "end", "strand",
 	 paste0(first, c("_coverage", "_methyl_coverage","_unmethyl_coverage")),
	 paste0(second, c("_coverage", "_methyl_coverage","_unmethyl_coverage")),
   	 "pvalue", "qvalue", "meth_diff"
	)
	
	all <- all[order(all\$qvalue), ]

	write.table(all, file=paste0(prefix, ".CpG.methylBase.summary.xls"),
		quote=F, sep="\t", row.names=F, col.names=T)

	qvalue <- 0.1
	difference <- 20

	all_diff <- all[all\$qvalue<=qvalue,]

	all_diff.hyper <- all_diff[all_diff\$meth_diff >= difference, ]
	all_diff.hypo <- all_diff[all_diff\$meth_diff <= (-1)*difference, ]

	write.table(all_diff.hyper, 
		file=paste0(prefix, ".CpG.methylBase.up.xls"),
		quote=F, sep="\t", row.names=F, col.names=T)

	write.table(all_diff.hypo, 
		file=paste0(prefix,".CpG.methylBase.dw.xls"),
		quote=F, sep="\t", row.names=F, col.names=T)


	png(paste0(prefix, ".CpG.bimodal.normal.distribution.png"))
	myMixmdl <- myDiff.to.mixmdl(myDiff,  plot=T,  main="C_9.vs.T_6")
	dev.off()

	png(paste0(prefix, ".CpG.cost.function.png"))
	plotCost(myMixmdl,  main="Cost function")
	dev.off()

	mydmr=edmr(myDiff, DMC.qvalue=0.1, DMC.methdiff=20, mode=1, ACF=TRUE)
	mysigdmr=filter.dmr(mydmr, DMR.qvalue=0.1, mean.meth.diff=20, 
		num.CpGs = 5, num.DMCs=3)

	write.table(as.data.frame(mysigdmr), file=paste0(prefix, ".CpG.eDMR.summary.xls"),
		quote=F, sep="\t", row.names=F, col.names=T)
		
	genebody=genebody.anno(file="/MPATHB/resource/UCSC/oar3.1/anno/oar3.1.all.type.bed")
	cpgi=cpgi.anno(file="/MPATHB/resource/UCSC/oar3.1/anno/oar3.1.CpGisland.bed")

	png(paste0(prefix, ".eDMR_genebody_annotation.png"))
	plotdmrdistr(mysigdmr,  genebody)
	dev.off()

	png(paste0(prefix, ".eDMR_cpgi_annotation.png"))
	plotdmrdistr(mysigdmr, cpgi)
	dev.off()

	dmr.genes=get.dmr.genes(myDMR=mysigdmr, subject=genebody\$promoter,id.type="gene.symbol")
	
	return(list(myDiff,mydmr,mysigdmr))
}

if (! ${called}) {
	print("Read in sam files and perform methylation calls")
	methylObj <- read.bismark(bamL, sample.id=name, assembly="${assembl}",
		save.context=typeL, read.context="none",
		mincov=${min_read_cov}, minqual=20,
		save.folder="${output_dir}/${min_read_cov}", 
		nolap=TRUE, treatment=group)
		for (type in typeL){
			resultList <- dmr(type, pbj, name)
		}
} else {
	print("Read in count files and perform methylation calls")
	pbj <- read(bamL, sample.id=name, assembly="${assembl}",
		treatment=group)
	resultList <- dmr('CpG', pbj, name)
}

#for (type in typeL) {
#	print(paste0("Analyzing ", type))
#	file_list <- paste0("${output_dir}/", name, "_", type, ".txt")
#
#	pbj <- read(file_list, sammple.id=name,assembly="${assembl}",
#		treatmeant=group)
#	if (type == "CpG") {
#		meth <- unite(pbj, destrand=TRUE)
#	} else {
#		meth <- unite(pbj, destrand=FALSE)
#	}
#	print("----Begin correlation plot")
#
#	print("Descriptive statistics")
#}

END


if [ "$execute" == "TRUE" ]; then
	Rscript --save ${output_dir}/${prefix}${mid}.r
#if [ "$?" == "0" ]; then /bin/rm -f ${file}${mid}.r; fi
fi

