#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:
    This is designed to generate Monocle Rmd file for doing scRNA-seq analysis.

JSON file:

{
	"ercc_comment": "定义有没有ERCC spike-in。如果有ERCC-spike-in，赋值为 TRUE, 则ERCC相关的代码块会执行。", 
	"ercc":"FALSE", 

	"UMI_comment": "定义有没有使用UMI。如果有UMI，赋值为TRUE，则UMI相关的代码块会执行。", 
	"UMI":"FALSE", 

	"species_comment":"定义物种",
	"species":"human",

	"nmads_comment": "筛选参数;数字越大，过滤掉的样品越少", 
	"nmads_seq_depth":3, 
	"nmads_gene_num":3,
	"nmads_mito_percent":3,
	"nmads_ercc_percent":3,


	"output_dir_comment": "数据目录 (ending slash needed)", 
	"output_dir": "/MPATHC/ct/ehbio/nankai_liulin_tumor/new/result/",

	"link_dir_comment": "不修改", 
	"link_dir": "result/", 

	"prefix_comment": "输出文件前缀", 
	"prefix":"tumor831", 


	"cnt_file": "/MPATHC/ct/ehbio/nankai_liulin_tumor/new/expr_count", 
	"pheno_file": "/MPATHC/ct/ehbio/nankai_liulin_tumor/new/phenoData", 
	"feature_file": "/MPATHC/ct/ehbio/nankai_liulin_tumor/new/featureData"

}

'''

import sys
import os
from json import dumps as json_dumps
from json import load as json_load
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
reload(sys)
sys.setdefaultencoding('utf8')

debug = 0

def fprint(content):
    """ 
    This is a Google style docs.

    Args:
        param1(str): this is the first param
        param2(int, optional): this is a second param
            
    Returns:
        bool: This is a description of what is returned
            
    Raises:
        KeyError: raises an exception))
    """
    print json_dumps(content,indent=1)

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print >>sys.stderr, desc
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--config-json", dest="config_json",
        help="A JSON file containging all input parameters.")
#    parser.add_option("-m", "--count-matrix-file", dest="matrix_file",
#        help="A matrix file according to format specified above. Absolute path should be used.")
#    parser.add_option("-p", "--pheno-file", dest="pheno_file",
#        help="(lowercase p) A pheno file according to format specified above. Absolute path should be used.")
#    parser.add_option("-f", "--feature-file", dest="feature_file",
#        help="A feature file according to format specified above. Absolute path should be used.")
#    parser.add_option("-o", "--output-dir", dest="output_directory",
#        help="Specify the output directory of result file.")
#    parser.add_option("-L", "--link-dir", dest="link_directory",
#        help="Specify the linking directory of result file.")
#    parser.add_option("-P", "--prefix", dest="prefix",
#        help="(uppercase P) Prefix for output files.")
#    parser.add_option("-e", "--ercc", dest="ercc",
#        help="Do the data contain ercc (TRUE) or not (FALSE).")
#    parser.add_option("-u", "--umi", dest="umi",
#        help="Is the data containing UMI (TRUE) or not (FALSE).")
    parser.add_option("-t", "--tag", dest="tag", default='monocle', 
        help="A string to discriminate this analysis with others. Default <monocle>.")
    parser.add_option("-d", "--datatype", dest="datatype",
        help="<UMI_count> (not the same count that given to DESeq2) or <FPKM_TPM> (Raw FPKM or TPM value)")
    parser.add_option("-f", "--filter", dest="filter_data",
        default="FALSE", 
        help="Default <FALSE>. If <TRUE>, filter data with default parameters.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.config_json != None, "A label needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    config_json = options.config_json
    config_jsonD = json_load(open(config_json))
    datatype = options.datatype
    tag = options.tag
    filter_data = options.filter_data
    config_jsonD['datatype'] = config_jsonD.get('datatype', datatype)
    config_jsonD['tag'] = config_jsonD.get('tag', tag)
    config_jsonD['filter_data'] = config_jsonD.get('filter_data', filter_data)
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    print '''# Monocle数据导入 {{#monocle-preprocess-{config_jsonD[tag]} }}

```{{r {config_jsonD[tag]}-setup2monocle}}
if(exists("debug")){{
  debug=FALSE
}} else {{
  debug=TRUE
}}
```


```{{r {config_jsonD[tag]}-installpackages,eval=F}}
ource("http://bioconductor.org/biocLite.R")
#biocLite()
biocLite(c('Biobase', 'VGAM', 'DDRTree', 'fastICA', 'matrixStats', 'densityClust', 'Rtsne', 'limma', 'qlcMatrix', 'proxy', 'slam' ))
#devtools::install_github("cole-trapnell-lab/monocle-release@develop")
biocLite("monocle")
```


```{{r {config_jsonD[tag]}-load-packages, include=FALSE, eval=TRUE}}
library(monocle)
library(reshape2)
options(stringsAsFactors = FALSE)
```

```{{asis echo=debug}}
## 数据输入文件和参数配置 {{#config}}
```

```{{r {config_jsonD[tag]}-parameter}}
# 需要修改的参数
#
# 定义有没有ERCC spike-in。如果有ERCC-spike-in，赋值为 TRUE, 则ERCC相关的代码块会执行。
{config_jsonD[tag]}_ercc <- {config_jsonD[ercc]}

# 定义有没有使用UMI。如果有UMI，赋值为TRUE，则UMI相关的代码块会执行。
{config_jsonD[tag]}_UMI <- {config_jsonD[UMI]}

# 定义物种
species <- '{config_jsonD[species]}'

# 筛选参数
# #数字越大，过滤掉的样品越少
{config_jsonD[tag]}_nmads_seq_depth={config_jsonD[nmads_seq_depth]}
{config_jsonD[tag]}_nmads_gene_num={config_jsonD[nmads_gene_num]}
{config_jsonD[tag]}_nmads_mito_percent={config_jsonD[nmads_mito_percent]}
{config_jsonD[tag]}_nmads_ercc_percent={config_jsonD[nmads_ercc_percent]}


# 数据目录 (ending slash needed)
{config_jsonD[tag]}_output_dir <- "{config_jsonD[output_dir]}"

# 数据用于报告的链接目录，相对目录，不含子目录
{config_jsonD[tag]}_link_dir <- "{config_jsonD[link_dir]}"

# 输出文件前缀
{config_jsonD[tag]}_prefix <- "{config_jsonD[prefix]}"

# 根据前述定义自动生成的参数和文件名字，一般不需要修改
# 需要注意的是，原始count数据文件默认存储在数据目录的summary目录，后缀为.rc.xls

# 数据文件
{config_jsonD[tag]}_cnt_file <- "{config_jsonD[cnt_file]}"
{config_jsonD[tag]}_pheno_file <- "{config_jsonD[pheno_file]}"
{config_jsonD[tag]}_feature_file <- "{config_jsonD[feature_file]}"

# 输出文件
{config_jsonD[tag]}_rds_file <- paste0({config_jsonD[tag]}_output_dir, {config_jsonD[tag]}_prefix, '.rds')
{config_jsonD[tag]}_rds_filter_file <- paste0({config_jsonD[tag]}_output_dir, {config_jsonD[tag]}_prefix, '.filter.rds')
{config_jsonD[tag]}_cnt_filter_file <- paste0({config_jsonD[tag]}_output_dir,{config_jsonD[tag]}_prefix,".filter_rc.xls")
{config_jsonD[tag]}_pheno_filter_file <- paste0({config_jsonD[tag]}_output_dir,{config_jsonD[tag]}_prefix,".filter_pheno.xls")

# 
batch_effect=FALSE
```

```{{r {config_jsonD[tag]}-mkdir-output-dir}}
system(paste0("mkdir -p ", {config_jsonD[tag]}_output_dir))
```

```{{r {config_jsonD[tag]}-biology-common-data}}
#定义线粒体基因
## MT-ATP6, MT-CYB, MT-ND1, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, MT-CO2, MT-CO1, MT-ND2, MT-ATP8, MT-CO3, MT-ND3
human_mtgene <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888", "ENSG00000198886", "ENSG00000212907", "ENSG00000198786", "ENSG00000198695", "ENSG00000198712", "ENSG00000198804", "ENSG00000198763","ENSG00000228253", "ENSG00000198938", "ENSG00000198840")

##
#mouse_mtgene <- c()

if(species == "human") {{
  mtgene <- human_mtgene
  cellcycle_rds <- "human_cycle_markers.rds"
}} else {{
  mtgene <- mouse_mtgene
  cellcycle_rds <- "mouse_cycle_markers.rds"  
}}
```

## 数据概览 {{#data_overview}}

```{{asis echo=debug}}

样品信息`phenoData`通常包含以下信息 (Table \@ref(tab:{config_jsonD[tag]}-pheno-example))。

Table: (\#tab:pheno-example) 样品信息数据，包括样品的名字、来源、类型、处理方式、提取时间、处理批次、上机批次等信息。信息的多少和类型可根据实际增删。这个表格的目的是对样品的关键信息进行注释，并且可以检测数据之间是否存在**批次效应**的影响。

|sample_id      |sample_type |individual |replicate | batch|well |seq_type | seq_lane| collection_time|
|:--------------|:-----------|:----------|:---------|-----:|:----|:--------|--------:|---------------:|
|WT.A.r1.1      |WT          |A          |r1        |     1|A01  |PE       |        1|        20161101|
|WT.B.r2.2      |WT          |B          |r2        |     2|A02  |PE       |        1|        20161102|
|WT.C.r3.3      |WT          |C          |r3        |     3|A03  |PE       |        1|        20161103|
|Induced.A.r1.1 |Induced.1h  |A          |r1        |     1|A04  |PE       |        1|        20161101|
|Induced.B.r2.2 |Induced.1h  |B          |r2        |     2|A05  |PE       |        1|        20161102|
|Induced.C.r3.3 |Induced.1h  |C          |r3        |     3|A06  |PE       |        1|        20161103|


基因信息数据`featureData`一般包括基因的名字 (ENSEMBLE ID, Gene symbol, ENtrez ID)、基因的注释、基因的长度、GC含量和在基因组上的位置，可以从[ENSEMBL BioMart](http://asia.ensembl.org/biomart/martview)下载。


单细胞基因表达数据以`SCESet`类的形式保存。这个类需要三个输入文件。

* exprs 表达量：表达值的数字矩阵，每一行为基因，每一列为样品。

* phenoData 样品信息, 每一行为样品名字，每一列为样品属性如细胞类型、培养条件、提取时间、处理方式等。

* featureData 基因信息, 每一行对应基因名字，每一列是基因的属性，如GeneSymbol、基因长度、GC含量、基因功能等。
```


```{{r {config_jsonD[tag]}-load-data, cache.extra=mtime(c(cnt_file, pheno_file, feature_file))}}
if(debug) {{
  print(cnt_file)
}}
{config_jsonD[tag]}_cnt_data <- read.table({config_jsonD[tag]}_cnt_file, sep="\\t", row.names=1, header=T,quote="", check.names=F)
{config_jsonD[tag]}_cnt_data <- {config_jsonD[tag]}_cnt_data[rowSums({config_jsonD[tag]}_cnt_data>0)>0,]

# Keep only samples appear in cnt_data and output in same order as in cnt_data.
{config_jsonD[tag]}_pheno_data <- read.table({config_jsonD[tag]}_pheno_file,sep="\\t", header=TRUE, row.names=1, quote="", check.names=F)

# Pay attention to drop=F to avoid aumotely transfeering single column dataframe to vector
{config_jsonD[tag]}_pheno_data <- {config_jsonD[tag]}_pheno_data[match(colnames({config_jsonD[tag]}_cnt_data), rownames({config_jsonD[tag]}_pheno_data)),, drop=F]
{config_jsonD[tag]}_pheno_data_col = colnames({config_jsonD[tag]}_pheno_data)
{config_jsonD[tag]}_pheno_data_col_len <- length({config_jsonD[tag]}_pheno_data_col)
if({config_jsonD[tag]}_pheno_data_col_len<2) {{
  {config_jsonD[tag]}_pheno_data_col <- rep({config_jsonD[tag]}_pheno_data_col,2)
}}

#if(({config_jsonD[tag]}_pheno_data_col[1]!="conditions") | ({config_jsonD[tag]}_pheno_data_col[2]!="batch")) {{
# stop("***Pleace check if phenoData has conditions and batch at the second and third column.***")
#}}


# Keep only features appeared in cnt_data and output in same order as in cnt_data.
{config_jsonD[tag]}_feature_data <- read.table({config_jsonD[tag]}_feature_file, sep="\\t", header=T, row.names=1,quote="", check.names=F)
{config_jsonD[tag]}_feature_data_rowname <- rownames({config_jsonD[tag]}_feature_data)

# Pay attention to drop=F to avoid aumotely transfeering single column dataframe to vector
{config_jsonD[tag]}_feature_data <- {config_jsonD[tag]}_feature_data[match(rownames({config_jsonD[tag]}_cnt_data), {config_jsonD[tag]}_feature_data_rowname),, drop=F]
{config_jsonD[tag]}_feature_data_col <- colnames({config_jsonD[tag]}_feature_data)

if(! "gene_short_name" %in% {config_jsonD[tag]}_feature_data_col){{
  stop("***Must have a column containing a list of unique gene symbols with column name as <gene_short_name>.***")
}}
```


基因表达数据{config_jsonD[datatype]}如 (Table \@ref(tab:{config_jsonD[tag]}-expr-cnt-part)) 所示。

```{{r {config_jsonD[tag]}-expr-cnt-part}}
knitr::kable(
  head({config_jsonD[tag]}_cnt_data[,1:7]), booktabs=TRUE,
  caption="基因表达Count值展示 (前6行，前7列)。"
)
```


样品信息、处理方式、测序信息等如 (Table \@ref(tab:{config_jsonD[tag]}-pheno-part))所示。

```{{r {config_jsonD[tag]}-pheno-part}}
knitr::kable(
  head({config_jsonD[tag]}_pheno_data), booktabs=TRUE,
  caption="样品信息、处理方式、测序信息等 (前6行)。"
)
```

```{{asis echo=debug}}
把基因表达数据、样品信息数据、基因注释数据整合生成`CellDataSet`数据集。

Although Monocle can be used with raw read counts,  these are not directly proportional to expression values unless you normalize them by length,  so some Monocle functions could produce nonsense results. If you don't have UMI counts,  We recommend you load up FPKM or TPM values instead of raw read counts.

`expressionFamily=negbinomial.size()` for UMI count data, transcript count. Raw count not normalized by transcript length may not be appraiate.

`expressionFamily=tobit()` for raw TPM and FPKM data. Do not do log transform.
```

```{{r}}
lowerDetectionLimitThresh <- 0.1
```

```{{r {config_jsonD[tag]}-monocle}}
{config_jsonD[tag]}_pheno_dataframe <- new("AnnotatedDataFrame", data={config_jsonD[tag]}_pheno_data)
{config_jsonD[tag]}_feature_dataframe <- new("AnnotatedDataFrame", data={config_jsonD[tag]}_feature_data)


if ("{config_jsonD[datatype]}" == "FPKM_TPM") {{
    {config_jsonD[tag]}_HSMM <- newCellDataSet(as.matrix({config_jsonD[tag]}_cnt_data), 
        phenoData = {config_jsonD[tag]}_pheno_dataframe, 
        featureData = {config_jsonD[tag]}_feature_dataframe, 
        lowerDetectionLimit=lowerDetectionLimitThresh, expressionFamily = tobit(Lower=lowerDetectionLimitThresh))
    
    # estimate RNA coutns
    {config_jsonD[tag]}_cnt_data <- relative2abs({config_jsonD[tag]}_HSMM, method="num_genes")   
}}

# make a new CellDataSet using the RNA counts
{config_jsonD[tag]}_HSMM <- newCellDataSet(as(as.matrix({config_jsonD[tag]}_cnt_data), "sparseMatrix"), 
    phenoData={config_jsonD[tag]}_pheno_dataframe, featureData = {config_jsonD[tag]}_feature_dataframe, 
    lowerDetectionLimit = lowerDetectionLimitThresh, expressionFamily = negbinomial.size())

# compute size factors for normalization
{config_jsonD[tag]}_HSMM <- estimateSizeFactors({config_jsonD[tag]}_HSMM)

# compute dispersion for differential expression analysis later
{config_jsonD[tag]}_HSMM <- estimateDispersions({config_jsonD[tag]}_HSMM)
```

移除在所有样品中表达值都为0的基因，获得用于下游分析的数据集如 (Table \@ref(tab:{config_jsonD[tag]}-remove-all-zero))所示。

```{{r {config_jsonD[tag]}-remove-all-zero}}
sta <- as.data.frame(dim({config_jsonD[tag]}_HSMM))
rownames(sta) <- c("Number of genes","Number of samples")
colnames(sta) <- NULL
knitr::kable(t(sta), booktabs=TRUE,
  caption="表达矩阵统计。"
)
```

列出部分基因的表达状态统计

```{{r {config_jsonD[tag]}-expr-sta-samp-cnt }}
{config_jsonD[tag]}_HSMM <- detectGenes({config_jsonD[tag]}_HSMM, min_expr = lowerDetectionLimitThresh)
knitr::kable(head(fData({config_jsonD[tag]}_HSMM)), booktabs=TRUE, caption="表达状态统计")
```

```{{r}}
expressed_genes <- row.names(subset(fData({config_jsonD[tag]}_HSMM), num_cells_expressed >= 3))
```

列出当前样品的信息


```{{r {config_jsonD[tag]}-pdata-2 }}
knitr::kable(head(pData({config_jsonD[tag]}_HSMM)), booktabs=TRUE, caption="样品信息")
```


```{{r, eval={config_jsonD[filter_data]} }}
pData({config_jsonD[tag]}_HSMM)$Total_mRNAs <- Matrix::colSums(exprs({config_jsonD[tag]}_HSMM))
{config_jsonD[tag]}_HSMM <- {config_jsonD[tag]}_HSMM[,pData({config_jsonD[tag]}_HSMM)$Total_mRNAs < 1e6]
{config_jsonD[tag]}_upper_bound <- 10^(mean(log10(pData({config_jsonD[tag]}_HSMM)$Total_mRNAs)) +
    2 * sd(log10(pData({config_jsonD[tag]}_HSMM)$Total_mRNAs)))
{config_jsonD[tag]}_lower_bound <- 10^(mean(log10(pData({config_jsonD[tag]}_HSMM)$Total_mRNAs)) -
    2 * sd(log10(pData({config_jsonD[tag]}_HSMM)$Total_mRNAs)))
    
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_gene_count_bound.pdf")
qplot(Total_mRNAs, data=pData({config_jsonD[tag]}_HSMM), geom="density") +
   geom_vline(xintercept={config_jsonD[tag]}_lower_bound) +
   geom_vline(xintercept={config_jsonD[tag]}_upper_bound)
dev.off()
```


```{{r}}
{config_jsonD[tag]}_HSMM <- {config_jsonD[tag]}_HSMM[,pData({config_jsonD[tag]}_HSMM)$Total_mRNAs > {config_jsonD[tag]}_lower_bound &
              pData({config_jsonD[tag]}_HSMM)$Total_mRNAs < {config_jsonD[tag]}_upper_bound &
              pData({config_jsonD[tag]}_HSMM)$num_genes_expressed >= 500]
{config_jsonD[tag]}_HSMM <- detectGenes({config_jsonD[tag]}_HSMM, min_expr = lowerDetectionLimitThresh)
```

检查过滤后的基因的表达是否仍服从对数正态分布 (Figure \@ref(fig:{config_jsonD[tag]}-expr-log-normal-distrib))。

```{{r {config_jsonD[tag]}-expr-log-normal-distrib }}
# Log-transform each value in the expression matrix.
L <- log(exprs({config_jsonD[tag]}_HSMM[expressed_genes,])+lowerDetectionLimitThresh)

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_Distribution_standard_gene_expression_distrib.pdf")
qplot(value, geom="density", data=melted_dens_df) +
   stat_function(fun = dnorm, size=0.5, color='red') +
   xlab("Standardized log(FPKM)") +
   ylab("Density")
dev.off()
```


```{{asis echo=debug}}

The first step is to decide which genes to use in clustering the cells. We could use all genes, but we'd be including a lot of genes that are not expressed at a high enough level to provide a meaningful signal. Including them would just add noise to the system. We can filter genes based on average expression level, and we can additionally select genes that are unusually variable across cells. These genes tend to be highly informative about cell state. 

The `setOrderingFilter` function marks genes that will be used for clustering in subsequent calls to `clusterCells`, although we will be able to provide other lists of genes if we want. The `plot_ordering_genes` function shows how variability (dispersion) in a gene's expression depends on the average expression across cells. The red line shows Monocle's expectation of the dispersion based on this relationship. The genes we marked for use in clustering are shown as black dots, while the others are shown as grey dots. 
```

```{{r}}
{config_jsonD[tag]}_disp_table <- dispersionTable({config_jsonD[tag]}_HSMM)
{config_jsonD[tag]}_unsup_clustering_genes <- subset({config_jsonD[tag]}_disp_table, mean_expression >= 0.2 & dispersion_empirical>=1*dispersion_fit)
{config_jsonD[tag]}_HSMM_variance <- setOrderingFilter({config_jsonD[tag]}_HSMM, {config_jsonD[tag]}_unsup_clustering_genes$gene_id)
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_ordering_gene.pdf")
plot_ordering_genes({config_jsonD[tag]}_HSMM_variance)
dev.off()
```


```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_pc_variance_explained.pdf")
plot_pc_variance_explained({config_jsonD[tag]}_HSMM_variance, return_all = F, fastpath=F, norm_method="log", maxit=1000) 
dev.off()
```


```{{r}}
#tSNE
#DDRTree
{config_jsonD[tag]}_HSMM_variance <- reduceDimension({config_jsonD[tag]}_HSMM_variance, max_components=2, num_dim = 8,
                        reduction_method = 'tSNE', verbose = T, perplexity=5)
{config_jsonD[tag]}_HSMM_variance <- clusterCells({config_jsonD[tag]}_HSMM_variance, num_clusters=5)
```


```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.color_by_clusters.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, color_by="as.factor(Cluster)")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.follicleDIA1.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, 1, 2, color="follicleDIA")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.total_mRNAs.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, 1, 2, color="Total_mRNAs")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.num_genes_expressed.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, 1, 2, color="num_genes_expressed")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.stage.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, 1, 2, color="Stage")
dev.off()
```

排除掉检测到的基因数`~num_genes_expressed`或其它干扰因素的影响后，再聚类。干扰因素可从上面的聚类图中看出，若某一类因素与分类结果极其吻合但无明显生物学意义即可视为干扰因素，多个因素用`~factor1+factor2`表示。

```{{r, eval=F}}
batch_effect=TRUE
```

```{{r, eval=batch_effect}}
{config_jsonD[tag]}_HSMM_variance <- reduceDimension({config_jsonD[tag]}_HSMM_variance, max_components=2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T, 
                        residualModelFormulaStr = "~num_genes_expressed",
                        perplexity=5)
{config_jsonD[tag]}_HSMM_variance <- clusterCells({config_jsonD[tag]}_HSMM_variance, num_clusters=5)
```

```{{r, eval=batch_effect}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.stage2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, 1, 2, color="Stage")
dev.off()
```

```{{r, eval=batch_effect}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.color_by_clusters2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, color_by="as.factor(Cluster)")
dev.off()
```


```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.follicleDIA2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, 1, 2, color="follicleDIA")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_clusters.oocyteDIA2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_variance, 1, 2, color="oocyteDIA")
dev.off()
```

绘制发育曲线。

```{{r}}
if (batch_effect) {{
    {config_jsonD[tag]}_HSMM_variance <- reduceDimension({config_jsonD[tag]}_HSMM_variance, max_components=2, 
                        reduction_method = 'DDRTree', verbose = T, 
                        residualModelFormulaStr = "~num_genes_expressed")
}} else {{
    {config_jsonD[tag]}_HSMM_variance <- reduceDimension({config_jsonD[tag]}_HSMM_variance, max_components=2, 
                        reduction_method = 'DDRTree', verbose = T)
}}
{config_jsonD[tag]}_HSMM_variance <- orderCells({config_jsonD[tag]}_HSMM_variance)
#{config_jsonD[tag]}_HSMM_variance <- orderCells({config_jsonD[tag]}_HSMM_variance, reverse=T)
```

```{{r, eval=root_tracjectory}}
{config_jsonD[tag]}_HSMM_variance_GM_state <- function(cds){{
  if (length(unique(pData(cds)$State)) > 1){{
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }}else {{
    return (1)
  }}
}}
```

```{{r, eval=root_tracjectory}}
{config_jsonD[tag]}_HSMM_variance <- orderCells({config_jsonD[tag]}_HSMM_variance, \
    root_state={config_jsonD[tag]}_HSMM_variance_GM_state({config_jsonD[tag]}_HSMM_variance))
```


```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_trajectory_colorState.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_variance, color_by="State")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_trajectory_colorPseudotime.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_variance, color_by="Pseudotime")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_trajectory_colorPseudotime_facetState.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_variance, color_by="Pseudotime") + facet_wrap(~State, nrow=1)
dev.off()
```


```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_trajectory_colorStage.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_variance, color_by="Stage")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_cell_trajectory_colorfollicleDIA.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_variance, color_by="follicleDIA")
dev.off()
```

基于PseudoTime的差异基因鉴定

```{{r}}
{config_jsonD[tag]}_HSMM_variance_diff_test_res <- differentialGeneTest({config_jsonD[tag]}_HSMM_variance,
    fullModelFormulaStr="~sm.ns(Pseudotime)", cores=40)
```

```{{r}}
{config_jsonD[tag]}_HSMM_variance_diff_test_res2 <- diff_test_res[order({config_jsonD[tag]}_HSMM_variance_diff_test_res$pval), ]
{config_jsonD[tag]}_HSMM_variance_diff_test_res2[,c("gene_short_name", "pval", "qval")]
```

```{{r}}
{config_jsonD[tag]}_HSMM_variance_gene_subset <- rownames(head({config_jsonD[tag]}_HSMM_variance_diff_test_res2))
```

```{{r}}
plot_gene_jitter({config_jsonD[tag]}_HSMM_variance[{config_jsonD[tag]}_HSMM_variance_gene_subset, ], grouping="State", min_expr=0.1, plot_trend=T)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_genes_jitter_pseudotime.pdf")
plot_genes_in_pseudotime({config_jsonD[tag]}_HSMM_variance[{config_jsonD[tag]}_HSMM_variance_gene_subset,], color_by="Stage")
dev.off()
```

筛选差异表达的基因

```{{r}}
{config_jsonD[tag]}_HSMM_variance_sig_gene_names <- row.names(subset({config_jsonD[tag]}_HSMM_variance_diff_test_res, qval < 0.01))
length(sig_gene_names)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_pseudotime_heatmap.pdf")
plot_pseudotime_heatmap({config_jsonD[tag]}_HSMM_variance[{config_jsonD[tag]}_HSMM_variance_sig_gene_names,],
                        num_clusters = 4,
                        cores = 10,
                        show_rownames = T)
dev.off()
```

Selecting genes based on PCA loading

```{{r}}
{config_jsonD[tag]}_HSMM_pca <- {config_jsonD[tag]}_HSMM[expressed_genes,]
{config_jsonD[tag]}_exprs_filtered <- t(t(exprs({config_jsonD[tag]}_HSMM_pca)/pData({config_jsonD[tag]}_HSMM_pca)$Size_Factor))
{config_jsonD[tag]}_exprs_filtered@x <- log({config_jsonD[tag]}_exprs_filtered@x+1)


# Calculate the variance across genes without converting to a dense matrix
{config_jsonD[tag]}_expression_means <- Matrix::rowMeans({config_jsonD[tag]}_exprs_filtered)
{config_jsonD[tag]}_expression_vars <- Matrix::rowMeans(({config_jsonD[tag]}_exprs_filtered-{config_jsonD[tag]}_expression_means)^2)

# Filter out genes that are constant across all cells
{config_jsonD[tag]}_genes_to_keep <- {config_jsonD[tag]}_expression_vars>0
{config_jsonD[tag]}_exprs_filtered <- {config_jsonD[tag]}_exprs_filtered[{config_jsonD[tag]}_genes_to_keep,]
{config_jsonD[tag]}_expression_means <- {config_jsonD[tag]}_expression_means[{config_jsonD[tag]}_genes_to_keep]

{config_jsonD[tag]}_expression_vars <- {config_jsonD[tag]}_expression_vars[{config_jsonD[tag]}_genes_to_keep]
```



```{{r}}

# Here's how to take the top PCA loading genes,  but using
# sparseMatrix operations the whole time,  using irlba. Note
# that the v matrix from irlba is the loading matrix

set.seed(0)

{config_jsonD[tag]}_irlba_pca_res <- irlba(t({config_jsonD[tag]}_exprs_filtered), 
    nu=0, center={config_jsonD[tag]}_expression_means, 
    scale=sqrt({config_jsonD[tag]}_expression_vars), right_only=T)$v

row.names({config_jsonD[tag]}_irlba_pca_res) <- row.names({config_jsonD[tag]}_exprs_filtered)

# Here, we will just
# take the top 200 genes from components 2 and 3.
# Component 1 usually is driven by technical noise.
# We could also use a more principled approach,
# similar to what dpFeature does below
{config_jsonD[tag]}_PC2_genes <- names(sort(abs({config_jsonD[tag]}_irlba_pca_res[, 2]), decreasing = T))[1:200]
{config_jsonD[tag]}_PC3_genes <- names(sort(abs({config_jsonD[tag]}_irlba_pca_res[, 3]), decreasing = T))[1:200]
{config_jsonD[tag]}_ordering_genes <- union({config_jsonD[tag]}_PC2_genes, {config_jsonD[tag]}_PC3_genes)
```

```{{r}}
{config_jsonD[tag]}_HSMM_pca <- setOrderingFilter({config_jsonD[tag]}_HSMM_pca, {config_jsonD[tag]}_ordering_genes)
#{config_jsonD[tag]}_HSMM_pca <- reduceDimension({config_jsonD[tag]}_HSMM_pca, max_components=2)
```

```{{r}}
#tSNE
#DDRTree
{config_jsonD[tag]}_HSMM_pca <- reduceDimension({config_jsonD[tag]}_HSMM_pca, max_components=2, num_dim = 8,
                        reduction_method = 'tSNE', verbose = T, perplexity=5)
{config_jsonD[tag]}_HSMM_pca <- clusterCells({config_jsonD[tag]}_HSMM_pca, num_clusters=5)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.color_by_clusters.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, color_by="as.factor(Cluster)")
dev.off()
```


```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.follicleDIA1.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, 1, 2, color="follicleDIA")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.total_mRNAs.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, 1, 2, color="Total_mRNAs")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.num_genes_expressed.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, 1, 2, color="num_genes_expressed")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.stage.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, 1, 2, color="Stage")
dev.off()
```

排除掉检测到的基因数`~num_genes_expressed`或其它干扰因素的影响后，再聚类。干扰因素可从上面的聚类图中看出，若某一类因素与分类结果极其吻合但无明显生物学意义即可视为干扰因素，多个因素用`~factor1+factor2`表示。

```{{r, eval=F}}
{config_jsonD[tag]}_HSMM_pca <- reduceDimension({config_jsonD[tag]}_HSMM_pca, max_components=2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T, 
                        residualModelFormulaStr = "~num_genes_expressed",
                        perplexity=5)
{config_jsonD[tag]}_HSMM_pca <- clusterCells({config_jsonD[tag]}_HSMM_pca, num_clusters=5)
```


```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.color_by_clusters.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, color_by="as.factor(Cluster)")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.stage2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, 1, 2, color="Stage")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.follicleDIA2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, 1, 2, color="follicleDIA")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_clusters.oocyteDIA2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_pca, 1, 2, color="oocyteDIA")
dev.off()
```

绘制发育曲线。

```{{r}}
if (batch_effect) {{
    {config_jsonD[tag]}_HSMM_pca <- reduceDimension({config_jsonD[tag]}_HSMM_pca, max_components=2, 
                        reduction_method = 'DDRTree', verbose = T, 
                        residualModelFormulaStr = "~num_genes_expressed")
}} else {{
    {config_jsonD[tag]}_HSMM_pca <- reduceDimension({config_jsonD[tag]}_HSMM_pca, max_components=2, 
                        reduction_method = 'DDRTree', verbose = T)
}}
{config_jsonD[tag]}_HSMM_pca <- orderCells({config_jsonD[tag]}_HSMM_pca)
#{config_jsonD[tag]}_HSMM_pca <- orderCells({config_jsonD[tag]}_HSMM_pca, root_state=GM_state({config_jsonD[tag]}_HSMM_pca))
#plot_cell_trajectory({config_jsonD[tag]}_HSMM_pca, color_by="Stage")
```

```{{r, eval=F}}
{config_jsonD[tag]}_HSMM_pca_GM_state <- function(cds){{
  if (length(unique(pData(cds)$State)) > 1){{
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }}else {{
    return (1)
  }}
}}
```

```{{r, eval=F}}
{config_jsonD[tag]}_HSMM_pca <- orderCells({config_jsonD[tag]}_HSMM_pca, \
    root_state={config_jsonD[tag]}_HSMM_pca_GM_state({config_jsonD[tag]}_HSMM_pca))
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_trajectory_colorCluster.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_pca, color_by="as.factor(Cluster)")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_trajectory_colorState.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_pca, color_by="State")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_trajectory_colorPseudotime.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_pca, color_by="Pseudotime")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_trajectory_colorPseudotime_facetState.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_pca, color_by="Pseudotime") + facet_wrap(~State, nrow=1)
dev.off()
```


```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_trajectory_colorStage.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_pca, color_by="Stage")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_cell_trajectory_colorfollicleDIA.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_pca, color_by="follicleDIA")
dev.off()
```

基于PseudoTime的差异基因鉴定

```{{r}}
{config_jsonD[tag]}_HSMM_pca_diff_test_res <- differentialGeneTest({config_jsonD[tag]}_HSMM_pca,
    fullModelFormulaStr="~sm.ns(Pseudotime)", cores=40)
```

```{{r}}
{config_jsonD[tag]}_HSMM_pca_diff_test_res2 <- diff_test_res[order({config_jsonD[tag]}_HSMM_pca_diff_test_res$pval), ]
{config_jsonD[tag]}_HSMM_pca_diff_test_res2[,c("gene_short_name", "pval", "qval")]
```

```{{r}}
{config_jsonD[tag]}_HSMM_pca_gene_subset <- rownames(head({config_jsonD[tag]}_HSMM_pca_diff_test_res2))
```

```{{r}}
plot_gene_jitter({config_jsonD[tag]}_HSMM_pca[{config_jsonD[tag]}_HSMM_pca_gene_subset, ], grouping="State", min_expr=0.1)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_genes_jitter_pseudotime.pdf")
plot_genes_in_pseudotime({config_jsonD[tag]}_HSMM_pca[{config_jsonD[tag]}_HSMM_pca_gene_subset,], color_by="Stage")
dev.off()
```

筛选差异表达的基因

```{{r}}
{config_jsonD[tag]}_HSMM_pca_sig_gene_names <- row.names(subset({config_jsonD[tag]}_HSMM_pca_diff_test_res, qval < 0.01))
length(sig_gene_names)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_pca_plot_pseudotime_heatmap.pdf")
plot_pseudotime_heatmap({config_jsonD[tag]}_HSMM_pca[{config_jsonD[tag]}_HSMM_pca_sig_gene_names,],
                        num_clusters = 4,
                        cores = 10,
                        show_rownames = T)
dev.off()
```


Selecting genes based on density peak clustering

Selecting genes for ordering cells using the techniques discussed above often works well in simple settings,  but for
more complex processes,  we recommend a new procedure called dpFeature. This procedure works by first projecting
cells into two dimensions using t-SNE and then detecting clusters using the\densityPeak" algorithm of Rodriguez and
Laio [4]. The resulting clusters of cells are compared using Monocle's differential expression analysis functions to find
genes that distinguish them. The top of these are then used to order cells.
To use dpFeature,  we first select superset of feature genes as genes expressed in at least 5% of all the cells.


```{{r}}
{config_jsonD[tag]}_HSMM_density <- detectGenes({config_jsonD[tag]}_HSMM, min_expr=0.1)
fData({config_jsonD[tag]}_HSMM_density)$use_for_ordering <- fData({config_jsonD[tag]}_HSMM_density)$num_cells_expressed>0.05*ncol({config_jsonD[tag]}_HSMM_density)
```

Then we will perform a PCA analysis to identify the variance explained by each PC (principal component). We can
look at a scree plot and determine how many pca dimensions are wanted based on whether or not there is a significant
gap between that component and the component after it. By selecting only the high loading PCs,  we ectively only
focus on the more interesting biological variations.

#look at the plot and decide how many di-
mensions you need. It is determined by a huge drop of variance at that dimension. pass that num-
ber to num_dim in the next function.

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_variance_plot_pc_variance_explained.pdf")
plot_pc_variance_explained({config_jsonD[tag]}_HSMM_variance, return_all = F, fastpath=F, norm_method="log", maxit=1000) 
dev.off()
```



```{{r}}
{config_jsonD[tag]}_HSMM_density <- reduceDimension({config_jsonD[tag]}_HSMM_density, 
    max_components=2, norm_method='log', num_dim=3, reduction_method='tSNE', verbose=T)
```

Then we can run density peak clustering to identify the clusters on the 2-D t-SNE space. Density peak algorithm
clusters cells based on each cell's local density () and the nearest distance config_jsonD of a cell to another cell with higher
distance. We can set a threshold for the ; config_jsonD and dene any cell with a higher local density and distance than the
thresholds as the density peaks. Those peaks are then used to dene the clusters for all cells. By default,  clusterCells
choose 95% of  and config_jsonD to dene the thresholds. We can also set a number of clusters (n) we want to cluster. In this
setting,  we will nd the top n cells with high config_jsonD with  among the top 50% range. The default setting often gives good
clustering.'

```{{r}}
{config_jsonD[tag]}_HSMM_density <- clusterCells({config_jsonD[tag]}_HSMM_density)
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.follicleDIA1.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color="follicleDIA")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.color_by_clusters.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, color_by="as.factor(Cluster)")
dev.off()
```


```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.total_mRNAs.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color="Total_mRNAs")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.num_genes_expressed.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color="num_genes_expressed")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.stage.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color="Stage")
dev.off()
```





排除掉检测到的基因数`~num_genes_expressed`或其它干扰因素的影响后，再聚类。干扰因素可从上面的聚类图中看出，若某一类因素与分类结果极其吻合但无明显生物学意义即可视为干扰因素，多个因素用`~factor1+factor2`表示。

```{{r, eval=F}}
{config_jsonD[tag]}_HSMM_density <- reduceDimension({config_jsonD[tag]}_HSMM_density, max_components=2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T, 
                        residualModelFormulaStr = "~num_genes_expressed",
                        perplexity=5)
{config_jsonD[tag]}_HSMM_density <- clusterCells({config_jsonD[tag]}_HSMM_density, num_clusters=5)
```


```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.stage2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color="Stage")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.cluster2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color_by="as.factor(Cluster)")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.follicleDIA2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color="follicleDIA")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.oocyteDIA2.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, 1, 2, color="oocyteDIA")
dev.off()
```

We also provide the decision plot for users to check the distribution of rho and delta to decide the threshold for defining the
cell clusters.

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.rho_delta.pdf")
plot_rho_delta({config_jsonD[tag]}_HSMM_density, rho_threshold=2, delta_threshold=4)
dev.off()
```

We could then re-run clustering based on the user defined threshold. To facilitate the computation,  we can set
(skip_rho_sigma = T) which enables us to skip the calculation of the rho and delta.

Selection of Rho and delta will determine the number of final clusters.

```{{r}}
{config_jsonD[tag]}_HSMM_density <- clusterCells({config_jsonD[tag]}_HSMM_density, 
        rho_threshold=2, delta_threshold=4, skip_rho_sigma=T, verbose=F)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_clusters.color_by_clusters_self_rho_delta.pdf")
plot_cell_clusters({config_jsonD[tag]}_HSMM_density, color_by="as.factor(Cluster)")
dev.off()
```

```{{r}}
{config_jsonD[tag]}_HSMM_density_clustering_DEG_genes <- differentialGeneTest({config_jsonD[tag]}_HSMM_density[expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 40)
```

```{{r}}
{config_jsonD[tag]}_HSMM_density_ordering_genes <- row.names({config_jsonD[tag]}_HSMM_density_clustering_DEG_genes)[order({config_jsonD[tag]}_HSMM_density_clustering_DEG_genes$qval)][1:1000]

{config_jsonD[tag]}_HSMM_density <- setOrderingFilter({config_jsonD[tag]}_HSMM_density, 
        ordering_genes={config_jsonD[tag]}_HSMM_density_ordering_genes)
```

```{{r}}
# reduction_method: DDRTree, ICA, tSNE, SimplePPT, L1-graph, SGL-tree
{config_jsonD[tag]}_HSMM_density <- reduceDimension({config_jsonD[tag]}_HSMM_density)
```

```{{r, eval=F}}
{config_jsonD[tag]}_HSMM_density <- reduceDimension({config_jsonD[tag]}_HSMM_density, max_components=2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T, 
                        residualModelFormulaStr = "~num_genes_expressed",
                        perplexity=5)
{config_jsonD[tag]}_HSMM_density <- clusterCells({config_jsonD[tag]}_HSMM_density, num_clusters=5)
```


绘制发育曲线。

```{{r}}
{config_jsonD[tag]}_HSMM_density <- orderCells({config_jsonD[tag]}_HSMM_density)
#{config_jsonD[tag]}_HSMM_density <- orderCells({config_jsonD[tag]}_HSMM_density, root_state=GM_state({config_jsonD[tag]}_HSMM_density))
#plot_cell_trajectory({config_jsonD[tag]}_HSMM_density, color_by="Stage")
```

```{{r, eval=F}}
{config_jsonD[tag]}_HSMM_density_GM_state <- function(cds){{
  if (length(unique(pData(cds)$State)) > 1){{
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }}else {{
    return (1)
  }}
}}
```

```{{r, eval=F}}
{config_jsonD[tag]}_HSMM_density <- orderCells({config_jsonD[tag]}_HSMM_density, \
    root_state={config_jsonD[tag]}_HSMM_density_GM_state({config_jsonD[tag]}_HSMM_density))
```


```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_trajectory_colorState.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_density, color_by="State")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_trajectory_colorPseudotime.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_density, color_by="Pseudotime")
dev.off()
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_trajectory_colorPseudotime_facetState.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_density, color_by="Pseudotime") + facet_wrap(~State, nrow=1)
dev.off()
```


```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_trajectory_colorStage.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_density, color_by="Stage")
dev.off()
```

```{{r, eval=F}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_cell_trajectory_colorfollicleDIA.pdf")
plot_cell_trajectory({config_jsonD[tag]}_HSMM_density, color_by="follicleDIA")
dev.off()
```

基于PseudoTime的差异基因鉴定

```{{r}}
{config_jsonD[tag]}_HSMM_density_diff_test_res <- differentialGeneTest({config_jsonD[tag]}_HSMM_density,
    fullModelFormulaStr="~sm.ns(Pseudotime)", cores=40)
```

```{{r}}
{config_jsonD[tag]}_HSMM_density_diff_test_res2 <- diff_test_res[order({config_jsonD[tag]}_HSMM_density_diff_test_res$pval), ]
{config_jsonD[tag]}_HSMM_density_diff_test_res2[,c("gene_short_name", "pval", "qval")]
```

```{{r}}
{config_jsonD[tag]}_HSMM_density_gene_subset <- rownames(head({config_jsonD[tag]}_HSMM_density_diff_test_res2))
```

```{{r}}
plot_gene_jitter({config_jsonD[tag]}_HSMM_density[{config_jsonD[tag]}_HSMM_density_gene_subset, ], grouping="State", min_expr=0.1)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_genes_jitter_pseudotime.pdf")
plot_genes_in_pseudotime({config_jsonD[tag]}_HSMM_density[{config_jsonD[tag]}_HSMM_density_gene_subset,], color_by="Stage")
dev.off()
```

筛选差异表达的基因

```{{r}}
{config_jsonD[tag]}_HSMM_density_sig_gene_names <- row.names(subset({config_jsonD[tag]}_HSMM_density_diff_test_res, qval < 0.01))
length(sig_gene_names)
```

```{{r}}
pdf("{config_jsonD[prefix]}_{config_jsonD[tag]}_density_plot_pseudotime_heatmap.pdf")
plot_pseudotime_heatmap({config_jsonD[tag]}_HSMM_density[{config_jsonD[tag]}_HSMM_density_sig_gene_names,],
                        num_clusters = 4,
                        cores = 10,
                        show_rownames = T)
dev.off()
```




```{{r}}
system(paste0("mkdir -p ", link_dir))
system(paste0("/bin/cp -up ", output_dir,"* ", link_dir))
```

'''.format(config_jsonD=config_jsonD)

    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        print >>sys.stderr,\
            "--Successful %s" % strftime(timeformat, localtime())

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print >>fh, "%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime)
    fh.close()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


