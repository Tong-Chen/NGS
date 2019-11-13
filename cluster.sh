#!/bin/bash
#############
#CT##########
#############

#set -x
set -e
#set -u
usage()
{
cat <<EOF
${txtcyn}
Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to do clustering using kmeans.

${txtred}Warning${txtrst}: This script can not deal well with one element cluster.
Please check it yourself.

${txtbld}OPTIONS${txtrst}:
	-f	Data file (with header line, the first column is the
 		rowname, tab seperated)${bldred}[NECESSARY]${txtrst}
		rowname must be unique. If not, use kmeans.uniq.py.
	-c	The number of cluster wanted.${bldred}[NECESSARY]${txtrst}
	-M	The cluster method to be used, choices are kmeans, pam (k-medoids cluster, default), clara (same algorithms as pam, but more fast for large data).
	-P	The way to preprocess data before clustering.
		[${txtred} 
		0		   : means no preprocess
		1          : means scale data
		2          : means using the difference value between current
		           : column and the one before
		3          : means values in each row are divided by the
		           : sqrt(sum(squares of values in this row))
	    4(default) : means first scale data and then using the difference 
				   : value between current column and the one before
		Both 2, 3, 4 are designed for clustering genes with same
		changing trend together.
		${txtrst}]
	-p	Select the optimum cluster number by elbow algorithm. 
		[${txtred}Default FALSE.${txtrst}]
		Accept TRUE. 
		When this is set to TRUE, a plot of within groups
		sum of squares is made. Acroread will open it, you choose the
		cluster number before which there is a sharp decline. Then
		close the opened pdf, you will be asked to input the chosed
		number. When this is TRUE, the number after -c is the expected
		maximum cluster number. All number of clusters ranges from 2
		to maximum cluster number will be calculated.
	-t	Give the maximum try number to run kmeans multiple times and
		choose the one with the least withinss.
		[${txtred}Default 10.${txtrst}]
	-S	Use scaled data for plot.
		[${txtred}Default FALSE.${txtrst}]
	-d	Get contour plot using density2d. 
		Only suitable for samples using numbers as names.
		[${txtred}Default FALSE.${txtrst}]
	-e	Evaluation cluster.[Default FALSE, accept TRUE]
	-x	The xlab of cluster.[default Value]
	-y	The ylab of cluster.[default Variable]
	-m	The title of the picture.[default '']
	-i	Install needed packages.[default FALSE, given TRUE once if there is
		error like <there is no package>]
	-l	Parameters for s-plot lines. 
EOF
}

file=
center=
preprocess=4
plotwithinss='FALSE'
try=10
evaluation='FALSE'
xlab='Value'
ylab='Variable'
mainT=''
ly=''
ist='FALSE'
method='pam'
scale_plot='FALSE'
density2d='FALSE'

while getopts "hf:c:d:P:p:t:e:x:y:m:M:S:l:i:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		c)
			center=$OPTARG
			;;
		P)
			preprocess=$OPTARG
			;;
		p)
			plotwithinss=$OPTARG
			;;
		t)
			try=$OPTARG
			;;
		e)
			evaluation=$OPTARG
			;;
		d)
			density2d=$OPTARG
			;;
		x)
			xlab=$OPTARG
			;;
		y)
			ylab=$OPTARG
			;;
		m)
			mainT=$OPTARG
			;;
		M)
			method=$OPTARG
			;;
		S)
			scale_plot=$OPTARG
			;;
		i)
			ist=$OPTARG
			;;
		l)
			ly=$OPTARG
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
mid=''

#if [ "$preprocess" == '1' ]; then
#	mid=${mid}'.scale'
#fi

#if [ "$preprocess" == '2' ]; then
#	mid=${mid}'.diff'
#fi


#if [ "$preprocess" == '3' ]; then
#	mid=${mid}'.trend'
#fi

#if [ "${scale_plot}" == "FALSE" ]; then
#	preprocess=1
#fi
#
if [ "$plotwithinss" == 'TRUE' ]; then
	cat <<EOF >${file}${mid}.$center.kmeans.chooseClusterNumber.r
data <- read.table(file="$file", sep='\t', header=T, row.names=1,
check.names=FALSE)
data <- as.matrix(data)
#Delate rows containing only zero
data <- data[rowSums(data==0)<ncol(data),]

preprocess=${preprocess}
if(ncol(data) <=2 && (preprocess==4 || preprocess==2)){
	preprocess=3
}

if(preprocess == 0){
	kdata <- data
}else if(preprocess == 1){
	kdata <- t(apply(data,1,scale))
}else if(preprocess == 2){
	kdata <- t(apply(data,1,diff))
}else if(preprocess == 3){
	norm_factors_for_each_row <- sqrt(apply(data^2, 1, sum))
	kdata <- data / norm_factors_for_each_row
}

wss <- (nrow(kdata)-1)*sum(apply(kdata,2,var))

for (i in 2:$center) wss[i] <- sum(kmeans(kdata, centers=i,
	iter.max=100, nstart=25)\$withinss)

pdf("${file}${mid}.$center.kmeans.chooseClusterNumber.pdf")
plot(1:$center, wss, type="b", xlab="Number of Clusters", ylab="Within
groups sum of squares")
dev.off()
	
EOF

Rscript --save ${file}${mid}.$center.kmeans.chooseClusterNumber.r
if [ $? == 0 ]; then
acroread ${file}${mid}.$center.kmeans.chooseClusterNumber.pdf

read -p ">>>Input the cluster number : " center
else
	echo "Wrong Rscript"
fi
fi

cat <<EOF >${file}${mid}.$center.${method}.r
if (${ist}){
	install.packages("cluster", repo="http://cran.us.r-project.org")
	install.packages("psych", repo="http://cran.us.r-project.org")
	install.packages("fpc", repo="http://cran.us.r-project.org")
}
library(cluster)
library(psych)
library(fpc)
data <- read.table(file="$file", sep='\t', header=T, row.names=1,
check.names=FALSE)
data <- as.matrix(data)
#Delete rows containing only zero
data <- data[rowSums(data==0)<ncol(data),]

preprocess=${preprocess}
if(ncol(data) <=2 && (preprocess==4 || preprocess==2)){
	preprocess=3
}

if(preprocess == 0){
	kdata <- data
}else if(preprocess == 1){
	kdata <- t(apply(data,1,scale))
}else if(preprocess == 2){
	kdata <- t(apply(data,1,diff))
}else if(preprocess == 3){
	norm_factors_for_each_row <- sqrt(apply(data^2, 1, sum))
	kdata <- data / norm_factors_for_each_row
}else if(preprocess == 4){
	kdata <- t(apply(data,1,scale))
	kdata <- t(apply(kdata,1,diff))
}


method = "${method}"
if (method == "pam" || method == "clara"){
	num_row = dim(data)[1] 
	if (num_row > 65536) {
		method = "clara"
	}
}

if ("${method}" == "kmeans"){
	print("Try ${method} for the first time.")
	fit <- ${method}(kdata, centers=$center, iter.max=100, nstart=25)
	withinss <- sum(fit\$withinss)
	print(paste("Get withinss for the first run", withinss))
	for (i in 1:$try) {
		tmpfit <- ${method}(kdata, centers=$center, iter.max=100, nstart=25)
		tmpwithinss <- sum(tmpfit\$withinss)
		print(paste(("The additional "), i, 'run, withinss', tmpwithinss))
		if (tmpwithinss < withinss){
			withins <- tmpwithinss
			fit <- tmpfit
		}
	}
	fit_cluster = fit\$cluster
} else if (method == "pam") {
	fit_cluster <- pam(kdata, k=$center, do.swap=FALSE, cluster.only=TRUE)
} else if (method == "clara") {
	num_row = dim(data)[1]	
	fit <- clara(kdata, k=$center, samples=num_row/$center)
	fit_cluster <- fit\$clustering
}

if ("${scale_plot}"){
	colnames_data = colnames(data)
	data <- t(apply(data,1,scale))
	colnames(data) = colnames_data
}

print("Output the mean value of cluster")
cluster.mean <- aggregate(data, by=list(fit_cluster), FUN=mean)
cluster.mean.colnames <- colnames(cluster.mean)
#cluster.mean.colnames[1] = paste('#',cluster.mean.colnames[1], sep='')
#colnames(cluster.mean) <- cluster.mean.colnames
write.table(t(cluster.mean), file="${file}${mid}.$center.${method}.cluster.mean.xls", sep='\t',col.names=F, row.names=T, quote=F)
print("Output the total sorted cluster name")
clust.out <- fit_cluster
kclust <- as.matrix(clust.out)
kclust.out <- cbind(kclust, data)
#means of n points in each cluster
mns <- sapply(split(data, fit_cluster), function(x) mean(unlist(x)))
#order the data
data.order <- data[order(order(mns)[fit_cluster]),]
write.table(data.order, file="${file}${mid}.$center.${method}.result", 
sep="\t", row.names=T, col.names=T, quote=F)

cluster <- clust.out
cluster <- as.data.frame(cluster)

dataWithClu <- cbind(ID=rownames(data), data, cluster)
dataWithClu <- dataWithClu[order(dataWithClu\$cluster),]
write.table(as.data.frame(dataWithClu), 
	file="${file}${mid}.${center}.${method}.result.final", 
	sep="\t", row.names=F, col.names=T, quote=F)

cluster_out <- cbind(ID=rownames(cluster), cluster)
cluster_out <- cluster_out[order(cluster_out\$cluster),]
write.table(as.data.frame(cluster_out), 
	file="${file}${mid}.${center}.${method}.result.factor_labeling", 
	sep="\t", row.names=F, col.names=F, quote=F)

#cluster_
#for(i in cluster){
#	ids <- 
#}

print("Plot the cluster data") 
png("${file}${mid}.$center.${method}.result.png")
xlabels <- colnames(data)
ylabels <- rev(rownames(data))
colN <- ncol(data.order)
rowN <- nrow(data.order)
image(1:colN, 1:rowN, t(data.order), xlab="$xlab",
ylab="$ylab", main="$mainT", xaxt='n', las=2)
axis(1, at=1:colN, labels=xlabels)
axis(2, at=1:rowN, labels=ylabels)
dev.off()
print("Output by cluster")
for(i in 1:$center){
	cluster <- data[kclust.out[,1]==i,]
	if(is.vector(cluster)){ #means only one value
		cluster <- t(as.matrix(cluster)) #have not tested
	}
	file_clu = paste("${file}${mid}.$center.${method}.cluster", i, sep='.')
	write.table(cluster, file=file_clu, sep="\t", row.names=T,
	col.names=T, quote=F)
}
if($evaluation){
print("Principle components plot")
png("${file}${mid}.$center.${method}.pc.png")
data.p <- prop.table(data, 1)*100
clusplot(data.p, fit_cluster, shade=F, labels=5, lines=0, color=T,
lty=4, main='Principal Components plot showing K-means clusters')
dev.off()
png("${file}${mid}.$center.${method}.png")
clusplot(data, fit_cluster, shade=T, labels=5, lines=0, color=T,
lty=4, main='K-means clusters')
dev.off()
png("${file}${mid}.$center.${method}.centroid.png")
plotcluster(data, fit_cluster)
dev.off()
kclust.out.p <- prop.table(as.matrix(kclust.out),1)*100
out <- capture.output(describe.by(kclust.out.p, kclust))
cat(out, file="${file}${mid}.$center.${method}.result.pc", sep='\n', append=F)
}
EOF

for i in `ls | grep "${file}${mid}.$center.${method}.cluster.[0-9]"`; do
	/bin/rm -f $i
done

Rscript ${file}${mid}.$center.${method}.r
echo 1>&2 "Please pay attention to one element cluster and correct
the result if any."
if [ $? == 0 ]; then
	s-plot lines -f ${file}${mid}.$center.${method}.cluster.mean.xls -B 0.5 -P top -o TRUE -R 45 -x "Samples" -c TRUE -C "rainbow($center)" -y "FPKM" 
dir=$(dirname ${file})
for i in `ls ${dir} | grep "${file}${mid}.$center.${method}.cluster.[0-9]"`; do
	awk 'BEGIN{OFS="\t"}{if(NR==1){print "Sample",$0}else print $0}' ${dir}/$i >${dir}/$i.tmp
	transpose.py ${dir}/$i.tmp >${dir}/$i.lines.xls
	s-plot lines -f ${dir}/$i.lines.xls -P none -B 0.5 -o TRUE -R 45 -x "Samples" -y "FPKM" 
	if [ "${density2d}" == "TRUE" ]; then
		s-plot density2d -f ${dir}/$i.lines.xls -o TRUE -y "Relative gene expression" 
	fi
	/bin/rm -f ${dir}/$i.tmp
done
else
	echo "Wrong"
fi

echo "Generate heatmap"

parseHeatmapSoutput.2.py -i \
	${file}${mid}.${center}.${method}.result.final \
	>${file}${mid}.${center}.${method}.result.final.sort.xls

s-plot heatmapS -f ${file}${mid}.${center}.${method}.result.final.sort.xls \
	-A 90 -T 1.5 -V 0.5 -l top -v 18 -u 15 -F 9 -j TRUE -M yellow -x green -y red -Z TRUE



/bin/mkdir -p ${file}${mid}.$center.${method}
/bin/mv -f ${file}${mid}.$center.${method}.* ${file}${mid}.$center.${method}
