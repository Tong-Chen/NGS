#!/bin/bash
#set -x
set -e
set -u
#########################################
#Data and format in <sample.cor>
#Control must befor treated sample.
#---here ends content(this line and '#' not included)---------------
#Cel file	Sample
#1.cel	control-1
#2.cel	control-2
#3.cel	control-3
#4.cel	MUT2-1
#5.cel	MUT2-2
#6.cel	MUT2-3
#---here ends content(this line and '#' not included)---------------
#Notes:The second column must not have same sample names
#########################################
#Notes:
#1.Chips based on different platform, like different GRL 
#number can not read together. 
#If this happens, error message 错误于read.affybatch(filenames = l$filenames, phenoData = l$phenoData,  : 
#  Cel file Name_UN-1.CEL does not seem to have the correct dimensions
#2.
#########################################

usage()
{
cat <<EOF
${txtcyn}
Usage (the least parameter):

$0 -f pheno.txt -s prefix${txtrst}

${bldblu}Function${txtrst}:

This script is used to analyze affy microarray.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, the first column is the
 		cel filename, the second column is the treat methods and
		different repeats (every string in this column should be
		unique), ${bldred}tab${txtrst} seperated)${bldred}[NECESSARY]
		[Control locates before treat]${txtrst}
		name    sample (this header line is needed)
		GSM723207.CEL   shCntrl1
		GSM723208.CEL   shCntrl2
		GSM723209.CEL   shCntrl3
		GSM723210.CEL   shDot1L1
		GSM723211.CEL   shDot1L2
		GSM723212.CEL   shDot1L3

	-s	The prefix of output file.
	-M	The algorithm to use.
		${bldred}Default [rma], accept mas5${txtrst}
	-i	If error happends when loading needed packages, plaease give
		TRUE to -i to install all needed ones.
		${bldred}Default [FALSE]${txtrst}
	-G	The file used to annotate affy probe ID.
		Currently accepted is a tow column file with the first column
		containing probe_sets and the second column containing gene
		symbols (header line is optional).
		${bldred}Not necessary but better to supply${txtrst}
	-m	The number of rows in layout. When drawing the raw
		image signal, this value is used.
	-n	The number of columns in layout. When drawing the raw 
		image signal, this value is used.
	-a	Execute PMA selection.[Default FALSE, accept TRUE]
	-r	Run the script[default] or only produce the script[FALSE].
	-g	Run DE test or not.[Default FALSE, accept TRUE only if there
		are two conditions in data. When TRUE,
		please supply other needed parameters.]
	-t	The repeat time of treated sample.[integer, default 3]
		{if -t and -c not setted, no PMA selection.
		This number relates with the number of rows in the last
		integer rows.
		}
	-c	The repeat time of control sample.[integer, default 3]
		{if -t and -c not setted, no PMA selection.
		This number relates with the number of rows in the first
		integer rows.
		}
	-p	The accepted maximum p-value.[default 0.05]
	-d	The accepted maximum fdr.[defaultault 0.3]
	-o	The accepted minimum fold change.
		[log2 based, default 1 means 2 times fold change.]
	-q	Need quality check, such as raw signals, pairwise comparision, 
		degradition, raw boxplot. Only needed when first run.
		[default:FALSE, skip this step. TRUE accepted if you run for 
		the first time.]
EOF
}

file=
prefix=
install='FALSE'
gpl=
pma='FALSE'
runDE='FALSE'
rawr=                #rows of raw signal
rawc=				 #cols of raw signal
run='TRUE'
controlR=3
treatR=3
pvalue=0.05
fdr=0.3
foldc=1
qcheck=FALSE
algorithm='rma'


while getopts "hf:s:i:M:m:n:a:r:g:G:t:c:p:d:o:q:" OPTION
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
			prefix=$OPTARG
			;;
		M)
			algorithm=$OPTARG
			;;
		i)
			install=$OPTARG
			;;
		m)
			rawr=$OPTARG
			;;
		n)
			rawc=$OPTARG
			;;
		a)
			pma=$OPTARG
			;;
		r)
			run=$OPTARG
			;;
		G)
			gpl=$OPTARG
			;;
		g)
			runDE=$OPTARG
			;;
		t)
			treatR=$OPTARG
			;;
		c)
			controlR=$OPTARG
			;;
		p)
			pvalue=$OPTARG
			;;
		d)
			fdr=$OPTARG
			;;
		o)
			foldc=$OPTARG
			;;
		q)
			qcheck=$OPTARG
			;;
		?)
			usage
			exit 1
			;;
	esac
done
if [ -z $file ] || [ -z $prefix ] ; then
	usage
	echo "No -f or -s supplied"
	exit 1
fi

if [ $qcheck == 'TRUE' ] && ( [ -z $rawr ] || [ -z $rawc ] ); then
	usage
	echo "No rawr or rawc for quality check."
	exit 1
fi

midname="affy."${algorithm}

cat <<EOF >$prefix.${midname}.r
if ($install){
	source("http://bioconductor.org/biocLite.R")
	#biocLite(c("affy", "genefilter", "gcrma", "affyPLM", 'xps', 'oligo'))
	biocLite(c("affy", "genefilter"))
	install.packages("amap", repo=" http://cran.us.r-project.org")
}

library("affy")
print("##reading $file")
pd <- read.AnnotatedDataFrame("$file", sep='\t', header=TRUE)
print('##reading cel files')
abc <- ReadAffy(filenames=rownames(pData(pd)), phenoData=pd, 
    sampleNames=pData(pd)[,1], verbose=TRUE)
if($qcheck){
	print('##Begin quality ckeck')
	#check the quality of microarray, the more uniform the better. 
	#or the less white points the better. Linux and windows may have 
	#different picture ways. Here windows preferred.
	#postscript will get a file very large, not recommended here.
	#postscript(file="$prefix.qc.raw.signal.eps", onefile=FALSE, 
	#	horizontal=FALSE, paper="special", width=10, height=8, pointsize=10)
	print('##Check raw signal')
	jpeg("$prefix.affy.qc.raw.signal.jpg", 2049, 2049, res=72, quality=95)  #width, height
	#the matrix 6 x 3, must change according to your sample
	par(mfrow=c($rawr,$rawc))  
	image(abc)
	dev.off() #needed to save the picture
	print('##Pairwise comparision')
	#Pairwise comparision. The more linear converge, the better.
	jpeg("$prefix.affy.MA.plot.jpg", 2048,2048,res=72,quality=95)
	MAplot(abc, pairs=TRUE, plot.method="smoothScatter")  #slow
	dev.off() #needed to save the picture
	#The boxplot, to see the expression variation and decide the normalization
	print('##Raw signal--boxplot')
	pdf("$prefix.affy.boxplot.raw.signal.pdf")
	#las is a parameter of par. 0: always parallel to the axis [_default_],
	#1: always horizontal, 2: always perpendicular to the axis, 
	#3: always vertical.
	#But when label is long, it will be cutted of.
	boxplot(abc, las=2, main="Boxplot of raw signal")  
	dev.off()
	print('##Begin RNA degradition plot')
	#QC-RNA degradation plots, mRNA degrades from 5' to 3'.
	deg <- AffyRNAdeg(abc)
	summaryAffyRNAdeg(deg)  #this is just a show of 'deg' 
	jpeg("$prefix.affy.rna.degradation.raw.signal.jpg", 800, 600,res=72,quality=95)
	plotAffyRNAdeg(deg, transform="shift.only")
	dev.off()
}
print("##Begin preprocessing data using << ${algorithm} >>")
#data preprocessing, background correction,normalization and summarization
#calculating expression
eset <- ${algorithm}(abc)

eset_expr <- exprs(eset)

if ("${algorithm}" == "mas5"){
	eset_expr <- log2(eset_expr)
}

#
print("##Normalized boxplot")
#The boxplot, to see the expression variation after normalization
pdf("$prefix.${midname}.boxplot.normalization.pdf")
boxplot(as.data.frame(eset_expr), las=2, main="Boxplot of raw signal")  
dev.off()
print("hclust of all samples")
#hcluster of RMA normalization
library("amap")
pdf("$prefix.${midname}.whole.gene.hcluster.pdf")
#t:returns the transpose of x.
plot(hcluster(t(eset_expr), method="pearson"), hang=-1)
dev.off()
print("#output expression value to $prefix.${midname}.expr")
#output
write.table(eset_expr, file="$prefix.${midname}.expr", sep="\t", 
	row.names=TRUE, col.names=TRUE, quote=FALSE)

if ($pma) {
	print("##Begin PMA detection")
	#pMA selection. if a probe is 'A' in all repeated sample no 
	#matter control or treated, delete this probe.
	detect <- exprs(mas5calls(abc, alpha1=0.01, alpha2=0.05))
	probe_A <- c()
	detectRC <- dim(detect)
	print('##annotate the probes which are A in all repeats')
	for(i in 1:detectRC[1]){
		delete <- 0
		treatDel <- 0
		half = $controlR
		for(j in 1:half){
			if(detect[i,j]=='A'){
				treatDel = treatDel + 1
			}
		}
		if(treatDel==3){
			delete <- 1
		}else{
			treatDel <- 0
			for(j in (half+1):(half+$treatR)){
				if(detect[i,j]=='A'){
					treatDel = treatDel + 1
				}
			}
			if(treatDel==3){
				delete <- 1
			}
		}
		if(delete==1){
			probe_A <- c(probe_A, 0) 
		}else{
			probe_A <- c(probe_A, 1)
		}
	}
	print('Output pma annotated result')
	detect <- cbind(detect, probe_A)  #only for testing
	write.table(detect, file="$prefix.${midname}.pma", sep="\t",
		row.names=TRUE, col.names=TRUE, quote=FALSE) 
	#combine the probe_A with espression data
	esetN <- cbind(eset_expr, probe_A)
	#delete untrustable data 
	filterdProbe <- subset(esetN, esetN[,'probe_A']==0)
	filterF <- filterdProbe[,1:detectRC[2]]
	esetNoA <- subset(esetN, esetN[,'probe_A']==1)
	#delete the last column
	esetF <- esetNoA[,1:detectRC[2]] #the final result for t-test
	print('output expression daate after filter and filtered probes')
	write.table(esetF, file="$prefix.${midname}.expr.aPMA", sep="\t", 
		row.names=TRUE, col.names=TRUE, quote=FALSE) 
	write.table(filterF, file="$prefix.${midname}.expr.fullA", sep="\t", 
		row.names=TRUE, col.names=TRUE, quote=FALSE) 
}else{
	esetF <- eset_expr
}

if ($runDE) {
	print('pick differentiated genes')
	library('genefilter')
	print('t-test')
	Ttest <- rowttests(esetF, as.factor(c(rep(1,$controlR),rep(2,$treatR))))
	print('method "BH" gives the false discovery rate ?p.adjust. We declare a collection of 100 genes with a maximum FDR of 0.10 to be differentially expressed (DE), then we expect a maximum of 10 genes to be false positives.')
	p.adjust <- p.adjust(Ttest\$p.value, method="BH")
	TtestAdj <- cbind(Ttest, p.adjust)
	esetFF <- cbind(esetF, TtestAdj)
	print("Output all gene expression with t-test results")
	write.table(esetFF,
		file="${prefix}${midname}_${foldc}_${pvalue}_${fdr}.expr.ttest",
		sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	#diffExpr <- subset(TtestAdjM, abs(TtestAdjM[,'dm'])>=$foldc &&
	#	TtestAdjM[,'p.adjust']<=$fdr && TtestAdjM[,'p.value']<=$pvalue)
	diffExpr <- subset(TtestAdj, abs(TtestAdj\$dm)>=$foldc)
	#diffExpr <- subset(esetFF,
	#abs(rowSums(esetFF[,2:$controlR])/rowSums(esetFF[,$controlR+1:$controlR+$treatR]))>=$foldc)
	diffExpr <- subset(diffExpr, diffExpr\$p.value<=$pvalue)
	diffExpr <- subset(diffExpr, diffExpr\$p.adjust<=$fdr)
	write.table(diffExpr, file="${prefix}${midname}_${foldc}_${pvalue}_${fdr}.deexpr",
		sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
}
EOF

if [ "$run" = 'TRUE' ];then
	Rscript $prefix.${midname}.r 
	/bin/rm -f $prefix.${midname}.r 
	sed -i '1 s/^/Probe\t/' ${prefix}.${midname}.expr
	#if [ -s ${gpl} ];then
	if ! [ -z ${gpl} ] && [ -s ${gpl} ];then
		#grep -v '^#' ${gpl} | cut -f 1,11,13 | \
		#	awk 'BEGIN{OFS="\t";FS="\t"}ARGIND==1{if(FNR>1)
		#	a[$1]=$1"@"$2"@"$3;}ARGIND==2{if(FNR>1) $1=a[$1]; print
		#	$0; }' - ${prefix}.expr | sed 's# /// #;#g' >${prefix}.expr.gene  
		#grep -v '^#' ${gpl} | paste - ${prefix}.expr >${prefix}.expr.anno
		awk 'BEGIN{OFS="\t";FS="\t"}ARGIND==1{
		a[$1]=$1"@"$2;}ARGIND==2{if(FNR>1) $1=a[$1]; print 
		$0; }' ${gpl} ${prefix}.${midname}.expr >${prefix}.${midname}.expr.gene  

	fi
fi


#cat <<EOF >$prefix.Makefile
#$prefix.r.CT:
#	@echo "The R script for analysis."
#
#$prefix.MA.plot.jpg.CT:
#	@echo "The pairwise comparision plot to indicate the
#	repeatability."
#$prefix.qc.raw.signal.jpg.CT:
#	@echo "The raw point picture to detect the quality of the
#	microarray. If no light points found, good. If run in linux,
#	light lines can be irgnored for system reason."
#
#$prefix.boxplot.raw.signal.pdf.CT:
#	@echo "The raw boxplot to see the difference among each sample
#	with raw data."
#
#$prefix.boxplot.normalization.pdf.CT:
#	@echo "The raw boxplot to see the difference among each sample
#	after normalization."
#
#$prefix.whole.gene.hcluster.pdf.CT:
#	@echo "The cluster of all samples. The tree graph can used
#	to tell if your result is right."
#
#$prefix.pma.CT:
#	@echo "The PMA detect result with 0 and 1 annotation in the last
#	column. 1 means kept. "
#
#$prefix.expr.CT:
#	@echo "The normalized expression value of all probes."
#
#$prefix.expr.fullA:
#	@echo "The expression value of wrong probes."
#
#$prefix.cor:
#	@echo "The input file of affy_full.sh."
#
#$prefix.t.test:
#	@echo "The filtered expression value with t-test results."
#
#${prefix}_${foldc}_${pvalue}_${fdr}.expr.ttest.CT:
#	@echo "The selected different expressed genes with expressed
#	values ttest results."
#
#${prefix}_${foldc}_${pvalue}_${fdr}.deexpr:
#	@echo "The selected different expressed genes and ttest results."
#EOF
