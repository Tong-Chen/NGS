#!/bin/bash

#set -x

usage()
{
cat <<EOF
${txtcyn}

***CREATED BY Chen Tong (chentong_biology@163.com)***

----Gene expr Matrix file--------------
# 常规表达矩阵，log2转换后或
## Deseq2的varianceStabilizingTransformation转换的数据
## 如果有批次效应，需要事先移除，可使用removeBatchEffect
## 如果有系统偏移(可用boxplot查看基因表达分布是否一致)，需要quantile normalization

Name	T0_1	T0_2	T0_3	T4_1	T4_2
TR19267|c0_g1|CYP703A2	1.431	0.77	1.309	1.247	0.485
TR19612|c1_g3|CYP707A1	0.72	0.161	0.301	2.457	2.794
TR60337|c4_g9|CYP707A1	0.056	0.09	0.038	7.643	15.379
TR19612|c0_g1|CYP707A3	2.011	0.689	1.29	0	0
TR35761|c0_g1|CYP707A4	1.946	1.575	1.892	1.019	0.999
TR58054|c0_g2|CYP707A4	12.338	10.016	9.387	0.782	0.563
TR14082|c7_g4|CYP707A4	10.505	8.709	7.212	4.395	6.103
TR60509|c0_g1|CYP707A7	3.527	3.348	2.128	3.257	2.338
TR26914|c0_g1|CYP710A1	1.899	1.54	0.998	0.255	0.427
----Matrix file--------------

----Trait information--------------

ID	WT	KO	OE Height Weight Diameter
T0_1	1	0	0	1	2	3
T0_2	1	0	0	2	4	6
T0_3	0	1	0	10	20	50
T4_1	0	1	0	15	30	80
T4_2	0	0	1	NA	9	8


Usage:

$0 options${txtrst}

${bldblu}Function${txtrst}:

This script is used to do WGCNA analysis and output desired .

The parameters for logical variable are either TRUE or FALSE.

${txtbld}OPTIONS${txtrst}:
	-f	Gene expression matrix (with header line, the first column is the
 		rowname, tab seperated. Colnames must be unique unless you
		know what you are doing. Normally rows should be genes and 
		columns should be samples.)${bldred}[NECESSARY]${txtrst}
	-t	Trait matrix file. Only numerical columsn will be used.
		${bldred}[Optional]${txtrst}
	-g	Specify the way to deal with duplicate row names.
		Default 0: representing duplicated row names are not allowed.
		Accept  1: representing make duplicated row names unique by adding <.1>, <.2> for the second, third appearances.
	-T	Network type. Default "unsigned", accept "signed" or "signed hybrid" (quotes needed).
	-c	Lowercase c. Correlation method used. Default "bicor", accept "pearson".
	-l	Log transform expression matrix.
		Default 0 (no transformation). Accept any non-zero positive number to
		initiate log2 transformation with given value as psuedo count to add to
		all numbers in gene expression matrix.
	-m	Merge Cut height. The larger value given, samller number of modules generated. [Default 0.2] 
	-e	Execute script (Default) or just output the script.
		[${bldred}Default TRUE${txtrst}]
	-i	Install the required packages. Normmaly should be TRUE if this is 
		your first time run s-plot.[${bldred}Default FALSE${txtrst}]
EOF
}

file=''
trait=''
mergeCutHeight=0.2
deal_duplicateRname=0
nettype="unsigned"
corType="bicor"
ist='FALSE'
log2=0
execute="TRUE"

while getopts "hf:t:g:T:c:l:e:i:" OPTION
do
	case $OPTION in
		h)
			echo "Help mesage"
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		t)
			trait=$OPTARG
			;;
		T)
			nettype=$OPTARG
			;;
		c)
			corType=$OPTARG
			;;
		l)
			log2=$OPTARG
			;;
		g)
			deal_duplicateRname=$OPTARG
			;;
		e)
			execute=$OPTARG
			;;
		i)
			ist=$OPTARG
			;;
		?)
			usage
			echo "Unknown parameters"
			exit 1
			;;
	esac
done

mid=".WGCNA"

if [ -z $file ] ; then
	echo 1>&2 "Please give filename."
	usage
	exit 1
fi



cat <<END >${file}${mid}.r

if ($ist){
	install.packages("pheatmap", repo="http://cran.us.r-project.org")
	source("https://bioconductor.org/biocLite.R")
	biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
	site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
	install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)
}

	
library(WGCNA)
library(reshape2)
library(stringr)

options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()

# 常规表达矩阵，log2转换后或
# Deseq2的varianceStabilizingTransformation转换的数据
# 如果有批次效应，需要事先移除，可使用removeBatchEffect
# 如果有系统偏移(可用boxplot查看基因表达分布是否一致)，
# 需要quantile normalization

exprMat <- "${file}"


# 官方推荐 "signed" 或 "signed hybrid"
# 用signed获得的模块包含的基因会少
type = "${nettype}"

# 相关性计算
# 官方推荐 biweight mid-correlation & bicor
# corType: pearson or bicor
corType = "${corType}"

corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

##导入数据##
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                     quote="", comment="", check.names=F)

dim(dataExpr)

head(dataExpr)[,1:8]

## ---- eval=T-------------------------------------------------------------
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                 max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg\$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg\$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg\$goodGenes], collapse = ",")));
  if (sum(!gsg\$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg\$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg\$goodSamples, gsg\$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)
head(dataExpr)[,1:8]

## ---- echo=T, fig.cap="查看是否有离群样品", fig.width=20-----------------
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
pdf(file="${file}.hcluster.tree.pdf", onefile=F, paper="special", 
	width=20, height=14, bg="white", pointsize=6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

## ---- echo=T-------------------------------------------------------------
## 软阈值筛选
## 软阈值的筛选原则是使构建的网络更符合无标度网络特征。
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)


pdf(file="${file}.softpower.pdf", onefile=F, paper="special", 
	width=22, height=14, bg="white", pointsize=6)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft\$fitIndices[,1], -sign(sft\$fitIndices[,3])*sft\$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft\$fitIndices[,1], -sign(sft\$fitIndices[,3])*sft\$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft\$fitIndices[,1], sft\$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft\$fitIndices[,1], sft\$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()

power = sft\$powerEstimate
power

## ---- echo=T-------------------------------------------------------------
# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
          ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
          ifelse(type == "unsigned", 6, 12))       
          )
          )
}

## ---- echo=T-------------------------------------------------------------
##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
# loadTOMs: 避免重复计算
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       networkType = type,
                       reassignThreshold = 0, mergeCutHeight = ${mergeCutHeight},
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOM=TRUE,
                       TOMDenom = "min",  deepSplit = 2,
                       detectCutHeight = 0.995, 
                       stabilityCriterion = "Individual fraction", 
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)


# 根据模块中基因数目的多少，降序排列，依次编号为 1-最大模块数。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net\$colors)

## ---- echo=T, fig.cap="层级聚类树展示各个模块"---------------------------
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net\$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
pdf(file="${file}.module_cluster.pdf", onefile=F, paper="special", 
	width=22, height=10, bg="white", pointsize=6)
plotDendroAndColors(net\$dendrograms[[1]], moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.5,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## ---- echo=T, fig.cap="模块之间的相关性", fig.width=5, fig.height=8------
### 基因和所在模块信息
gene_module <- data.frame(ID=colnames(dataExpr), module=moduleColors)
gene_module = gene_module[order(gene_module\$module),]
write.table(gene_module,file=paste0(exprMat,".gene_module.xls"),
            sep="\t",quote=F,row.names=F)

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net\$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

## 保存模块代表性信息
MEs_colt = as.data.frame(t(MEs_col))
colnames(MEs_colt) = rownames(dataExpr)
write.table(MEs_colt,file=paste0(exprMat,".module_eipgengene.xls"),
            sep="\t",quote=F)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距

pdf(file="${file}.moduleCorrelation.pdf", onefile=F, paper="special", 
	width=14, height=24, bg="white", pointsize=6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

## 如果有表型数据，也可以跟ME数据放一起，一起出图
#MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
#plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
#                      marDendro = c(3,3,2,4),
#                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
#                      xLabelsAngle = 90)


## ---- echo=T-------------------------------------------------------------
# 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# 否则需要再计算一遍，比较耗费时间
# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net\$TOMFiles[1], verbose=T)

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^power
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

## ---- echo=T, eval=F, fig.cap="TOM plot"---------------------------------
## # 这一部分特别耗时，行列同时做层级聚类
## TOMplot(plotTOM, net\$dendrograms, moduleColors,
##         main = "Network heatmap plot, all genes")

## ---- fig.cap="TOMplot"--------------------------------------------------
#knitr::include_graphics("images/TOMplot.png")

## ---- echo=T-------------------------------------------------------------
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整
cyt = exportNetworkToCytoscape(TOM,
             edgeFile = paste(exprMat, ".edges.txt", sep=""),
             nodeFile = paste(exprMat, ".nodes.txt", sep=""),
             weighted = TRUE, threshold = 0.01,
             nodeNames = probes, nodeAttr = moduleColors)

## ---- echo=T-------------------------------------------------------------
trait <- "${trait}"
# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]

	## ---- echo=T, fig.height=8, fig,width=8, fig.cap="模块与表型的关联"------
	### 模块与表型数据关联
	if (corType=="pearsoon") {
	  modTraitCor = cor(MEs_col, traitData, use = "p")
	  modTraitP = corPvalueStudent(modTraitCor, nSamples)
	} else {
	  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
	  modTraitCor = modTraitCorP\$bicor
	  modTraitP   = modTraitCorP\$p
	}
	# signif表示保留几位小数
	textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
	dim(textMatrix) = dim(modTraitCor)

	pdf(file="${file}.module_trait_correlation.pdf", onefile=F, paper="special", 
		width=22, height=20, bg="white", pointsize=6)
	labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
				   yLabels = colnames(MEs_col), 
				   cex.lab = 0.5, 
				   ySymbols = colnames(MEs_col), colorLabels = FALSE, 
				   colors = blueWhiteRed(50), 
				   textMatrix = textMatrix, setStdMargins = FALSE, 
				   cex.text = 0.5, zlim = c(-1,1),
				   main = paste("Module-trait relationships"))
	dev.off()

	modTraitCorMelt = as.data.frame(modTraitCor)
	write.table(modTraitCorMelt,file=paste0(exprMat,".module_trait_correlation.xls"),
				sep="\t",quote=F)
	modTraitCorMelt\$ID = rownames(modTraitCor)
	modTraitCorMelt = melt(modTraitCorMelt)
	colnames(modTraitCorMelt) <- c("Module","Trait","PersonCorrelationValue")
	modTraitPMelt = as.data.frame(modTraitP)
	write.table(modTraitPMelt,file=paste0(exprMat,".module_trait_correlationPvalue.xls"),
				sep="\t",quote=F)
	modTraitPMelt\$ID = rownames(modTraitP)
	modTraitPMelt = melt(modTraitPMelt)
	colnames(modTraitPMelt) <- c("Module","Trait","Pvalue")
	#modTraitCorP = cbind(modTraitCorMelt, Pvalue=modTraitPMelt\$Pvalue)
	modTraitCorP = merge(modTraitCorMelt, modTraitPMelt, by=c("Module","Trait"))
	write.table(modTraitCorP,file=paste0(exprMat,".module_trait_correlationPvalueMelt.xls"),
				sep="\t",quote=F,row.names=F)
}

## ---- echo=T, eval=F,  fig.width=8, fig.height=8-------------------------
## ## 从上图可以看到MEmagenta与Insulin_ug_l相关
## 
## ## 模块内基因与表型数据关联
## 
## # 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
## # 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
## # 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达
## # 值算出相关系数。
## # 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要
## # 。
## 
## ### 计算模块与基因的相关性矩阵
## 
## if (corType=="pearsoon") {
##   geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
##   MMPvalue = as.data.frame(corPvalueStudent(
##              as.matrix(geneModuleMembership), nSamples))
## } else {
##   geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
##   geneModuleMembership = geneModuleMembershipA\$bicor
##   MMPvalue   = geneModuleMembershipA\$p
## }
## 
## 
## # 计算性状与基因的相关性矩阵
## 
## ## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
## 
## if (corType=="pearsoon") {
##   geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
##   geneTraitP = as.data.frame(corPvalueStudent(
##              as.matrix(geneTraitCor), nSamples))
## } else {
##   geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
##   geneTraitCor = as.data.frame(geneTraitCorA\$bicor)
##   geneTraitP   = as.data.frame(geneTraitCorA\$p)
## }
## 
## 
## geneTraitCorMelt = as.data.frame(geneTraitCor)
## write.table(geneTraitCorMelt,file=paste0(exprMat,".gene_trait_correlation.xls"),
##             sep="\t",quote=F)
## geneTraitCorMelt\$ID = rownames(geneTraitCor)
## geneTraitCorMelt = melt(geneTraitCorMelt)
## colnames(geneTraitCorMelt) <- c("Gene","Trait","PersonCorrelationValue")
## geneTraitPMelt = as.data.frame(geneTraitP)
## write.table(geneTraitPMelt,file=paste0(exprMat,".gene_trait_correlationPvalue.xls"),
##             sep="\t",quote=F)
## geneTraitPMelt\$ID = rownames(geneTraitP)
## geneTraitPMelt = melt(geneTraitPMelt)
## colnames(geneTraitPMelt) <- c("Gene","Trait","Pvalue")
## #geneTraitCorP = cbind(geneTraitCorMelt, Pvalue=geneTraitPMelt\$Pvalue)
## geneTraitCorP = merge(geneTraitCorMelt, geneTraitPMelt, by=c("Gene","Trait"))
## write.table(geneTraitCorP,file=paste0(exprMat,".gene_trait_correlationPvalueMelt.xls"),
##             sep="\t",quote=F,row.names=F)
## 
## 
## # 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
## module = "magenta"
## pheno = "Insulin_ug_l"
## modNames = substring(colnames(MEs_col), 3)
## # 获取关注的列
## module_column = match(module, modNames)
## pheno_column = match(pheno,colnames(traitData))
## # 获取模块内的基因
## moduleGenes = moduleColors == module
## 
## sizeGrWindow(7, 7)
## par(mfrow = c(1,1))
## # 与性状高度相关的基因，也是与性状相关的模型的关键基因
## verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
##                    abs(geneTraitCor[moduleGenes, pheno_column]),
##                    xlab = paste("Module Membership in", module, "module"),
##                    ylab = paste("Gene significance for", pheno),
##                    main = paste("Module membership vs. gene significance\n"),
##                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

## ---- echo=T, eval=F-----------------------------------------------------
## ### 计算邻接矩阵
## adjacency = adjacency(dataExpr, power = power)
## 
## ### 把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵。
## TOM = TOMsimilarity(adjacency)
## dissTOM = 1-TOM
## 
## ### 层级聚类计算基因之间的距离树
## geneTree = hclust(as.dist(dissTOM), method = "average")
## 
## ### 模块合并
## # We like large modules, so we set the minimum module size relatively high:
## minModuleSize = 30
## # Module identification using dynamic tree cut:
## dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, cutHeight = 0.9,
##                             deepSplit = 2, pamRespectsDendro = FALSE,
##                             minClusterSize = minModuleSize)
## # Convert numeric lables into colors
## dynamicColors = labels2colors(dynamicMods)
## 
## ### 通过计算模块的代表性模式和模块之间的定量相似性评估，合并表达图谱相似
## #的模块
## MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
## MEs = MEList\$eigengenes
## # Calculate dissimilarity of module eigengenes
## MEDiss = 1-cor(MEs)
## # Cluster module eigengenes
## METree = hclust(as.dist(MEDiss), method = "average")
## MEDissThres = 0.25
## 
## # Call an automatic merging function
## merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
## # The merged module colors
## mergedColors = merge\$colors;
## # Eigengenes of the new merged
## 
## ## 分步法完结


END

if [ "$execute" == "TRUE" ]; then
	Rscript ${file}${mid}.r
	if [ "$?" == "0" ]; then 
		#/bin/rm -f ${file}${mid}.r
		/bin/rm -f Rplots.pdf	
	fi
fi

if test "${preprocess}" == "TRUE"; then
	/bin/mv -f ${file}".nostd0" ${file}
fi

#convert -density 200 -flatten ${file}${mid}.eps ${first}${mid}.png
