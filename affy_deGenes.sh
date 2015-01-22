#!/bin/bash
#set -x
set -e
set -u
#########################################

usage()
{
cat <<EOF
${txtcyn}
Usage (the least parameter):

$0 -f gene_expression.matrix -s prefix${txtrst}

${bldblu}Function${txtrst}:

This script is used to generate differentially expressed genes (DE
genes) using t-test or wilcox.test.

${txtbld}OPTIONS${txtrst}:
	-f	Data matrix file (with header line)
		${bldred}NECESSARY${txtrst}
	-R	The replication time for each group of samples.
		${bldred}NECESSARY, in format like '3,3,3' or '1,2,3,4'.
		The order matters.${txtrst}
	-v	The sample name of each group.
		${bldred}NECESSARY, in format like "'samp1','samp2','samp3'".
		Thr order matters.${txtrst}
	-s	Statistical method one wmat to used.
		${bldred}Defaulr t.test, accept wilcox.test.${txtrst}
	-p	The accepted maximum p-value.[default 0.05]
	-d	The accepted maximum fdr.[defaultault 0.3]
	-o	The accepted minimum fold change.
		[log2 based, default 1 means 2 times fold change.]
	-l	The data is already log2 transferred.
		${bldred}[Default TRUE, set to FALSE if not log2 transferred]${txtrst}
	-i	If error happends when loading needed packages, plaease give
		TRUE to -i to install all needed ones.
		${bldred}Default [FALSE]${txtrst}
	-r	Run the script[default] or only produce the script[FALSE].
EOF
}

file=
install='FALSE'
run='TRUE'
replication=
sampleName=
pvalue=0.05
fdr=0.3
foldc=1
log2='TRUE'
sta_m='t.test'


while getopts "hf:R:s:v:r:l:i:p:d:o:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		R)
			replication=$OPTARG
			;;
		v)
			sampleName=$OPTARG
			;;
		s)
			sta_m=$OPTARG
			;;
		l)
			log2=$OPTARG
			;;
		i)
			install=$OPTARG
			;;
		r)
			run=$OPTARG
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
		?)
			usage
			exit 1
			;;
	esac
done
if [ -z $file ] ; then
	usage
	echo "No -f supplied"
	exit 1
fi


midname="DE.${sta_m}"

cat <<EOF >$file.${midname}.r

library('genefilter')


run_DE <- function 
(a, 
 lena,
 data, 
 sampleName
){
	for(i in 1:(lena-1)){
		if (i>1) {
			start1 <- sum(a[1:(i-1)]) + 1
		} else {
			start1 <- 1
		}
		end1   <- start1 + a[i] - 1
		v1 <- c(rep(1, a[i]))
		samp1 <- sampleName[i]
		for(j in (i+1):lena){
			start2 <- sum(a[1:(j-1)]) +1
			end2   <- start2 + a[j] - 1
			end = sum(a[1:j])
			v2 <- c(rep(2, a[j]))
			#print(paste(start1, ":", end1, ";", start2, ":", end2,
			#			";"))
			#print(v1)
			#print(v2)
			samp2 <- sampleName[j]
			currentSample <- as.matrix(data[, c(start1:end1, start2:end2)])
			de_compute(currentSample, samp1, samp2, v1, v2)		
		}
	}
}


rowwilcox.test <- function(x, controlR, treatR) {
	p_fc_list <- c()
	inner_wilcox.test <- function (x, p_fc_list, controlR,treatR){
		control <- x[1:controlR]
		treat <- x[(controlR+1):(controlR+treatR)]
		p <- wilcox.test(control, treat)\$p.value
		fc <- mean(control) - mean(treat)
		if(length(p_fc_list)==0){
			p_fc_list <<- c(fc, p)
		}else {
			p_fc_list <<- rbind(p_fc_list, c(fc, p))
		}
	}
	apply(x, 1, function(x) inner_wilcox.test(x,p_fc_list,controlR,treatR))
	rownames(p_fc_list) <- rownames(x)
	return (as.data.frame(p_fc_list))
}

my_wilcox.test <- function
(esetF,
controlR, 
treatR
) {
	sta_test <- rowwilcox.test(esetF, controlR, treatR)
	colnames(sta_test) <- c('log2FC', 'p.value')
	return (sta_test)
}



de_compute <- function
(
esetF,
samp1, 
samp2,
v1,
v2
 ){
 	print(paste("Perform ${sta_m} for", samp1, "and", samp2))
	if ("${sta_m}" == "t.test"){
		Ttest <- rowttests(esetF, as.factor(c(v1,v2)))
		Ttest <- Ttest[, 2-3]
		colnames(Ttest) <- c('log2FC', 'p.value')
	} else if ("${sta_m}" == "wilcox.test"){
		Ttest <- my_wilcox.test(esetF, length(v1), length(v2))
	}
	p.adjust <- p.adjust(Ttest\$p.value, method="BH")
	TtestAdj <- cbind(Ttest, p.adjust)
	esetFF <- cbind(esetF, TtestAdj)
	sampC <- paste(samp1, samp2, sep="_")
	fileO <- paste("${file}","${midname}",sampC,"${foldc}","${pvalue}","${fdr}","expr", sep=".")
	write.table(esetFF, file=fileO, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	system(paste("sed -i '1 s/^/Gene\t/'", fileO))
	diffExpr <- subset(esetFF,  TtestAdj\$log2FC>=$foldc)
	diffExpr <- subset(diffExpr, diffExpr\$p.value<=$pvalue)
	diffExpr <- subset(diffExpr, diffExpr\$p.adjust<=$fdr)
	fileO <- paste("${file}","${midname}",sampC,"${foldc}","${pvalue}","${fdr}","deexpr.up", sep=".")
	write.table(diffExpr, file=fileO, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	system(paste("sed -i '1 s/^/Gene\t/'", fileO))
	foldc <- (-1) * ${foldc}
	diffExpr <- subset(esetFF,  TtestAdj\$log2FC<=foldc)
	diffExpr <- subset(diffExpr, diffExpr\$p.value<=$pvalue)
	diffExpr <- subset(diffExpr, diffExpr\$p.adjust<=$fdr)
	fileO <- paste("${file}","${midname}",sampC,"${foldc}","${pvalue}","${fdr}","deexpr.dw", sep=".")
	write.table(diffExpr, file=fileO, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	system(paste("sed -i '1 s/^/Gene\t/'", fileO))
}

data <- as.matrix(read.table(file="${file}", header=T, sep="\t",row.names=1))

if(! ${log2}){
	data <- log2(data)
}

replication <- c(${replication})
sampleName  <- c(${sampleName})

len_replication <- length(replication)

run_DE(replication, len_replication, data, sampleName)

EOF

if [ "$run" = 'TRUE' ];then
	Rscript $file.${midname}.r 
	/bin/rm -f $file.${midname}.r 
	#sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.expr.ttest
	#sed -i '1 s/^/Gene\t/' ${file}${midname}_${foldc}_${pvalue}_${fdr}.deexpr
fi


