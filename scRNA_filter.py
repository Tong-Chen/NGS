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
    This is designed to generate Rmd file for doing scRNA-seq filter.
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
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    print '''# 单细胞测序样品质控 {{#singlecell-preprocess}}

```{{r setup_2}}
if(exists("debug")){{
  debug=FALSE
}} else {{
  debug=TRUE
}}
```


```{{r install_packages, eval=F}}
#source('https://bioconductor.org/biocLite.R')
#biocLite('BiocInstaller')
#devtools::install_github(c("hemberg-lab/scRNA.seq.funcs","JustinaZ/pcaReduce","tallulandrews/M3D"))
#devtools::install_github("hms-dbmi/scde", build_vignettes = FALSE)
#install.packages(c("mvoutlier","ROCR"))
#biocLite(c("RUVSeq","pcaMethods","SC3","TSCAN","monocle","MultiAssayExperiment","SummarizedExperiment"))
#devtools::install_github("satijalab/seurat")
```


```{{r load_packages, include=FALSE, eval=TRUE}}
#source('https://bioconductor.org/biocLite.R')
#biocLite('BiocInstaller')
#devtools::install_github(c("hemberg-lab/scRNA.seq.funcs","JustinaZ/pcaReduce","tallulandrews/M3D"))
#devtools::install_github("hms-dbmi/scde", build_vignettes = FALSE)
#install.packages(c("mvoutlier","ROCR"))
#biocLite(c("RUVSeq","pcaMethods","SC3","TSCAN","monocle","MultiAssayExperiment","SummarizedExperiment"))
library(scater, quietly = TRUE)
library(Seurat,quietly=TRUE)
library(mvoutlier, quietly=TRUE)
#library(knitr, quietly = TRUE)
#library(scran, quietly = TRUE)
library(plyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library("gridExtra", quietly=TRUE)
library(sROC, quietly = TRUE)
#library("scde", quietly = TRUE)
#library("SC3", quietly=TRUE)
options(stringsAsFactors = FALSE)
```

```{{asis echo=debug}}
## 数据输入文件和参数配置 {{#config}}
```

```{{r parameter}}
# 需要修改的参数
#
# 定义有没有ERCC spike-in。如果有ERCC-spike-in，赋值为 TRUE, 则ERCC相关的代码块会执行。
ercc <- {config_jsonD[ercc]}

# 定义有没有使用UMI。如果有UMI，赋值为TRUE，则UMI相关的代码块会执行。
UMI <- {config_jsonD[UMI]}

# 定义物种
species <- '{config_jsonD[species]}'

# 筛选参数
# #数字越大，过滤掉的样品越少
nmads_seq_depth={config_jsonD[nmads_seq_depth]}
nmads_gene_num={config_jsonD[nmads_gene_num]}
nmads_mito_percent={config_jsonD[nmads_mito_percent]}
nmads_ercc_percent={config_jsonD[nmads_ercc_percent]}


# 数据目录 (ending slash needed)
output_dir <- "{config_jsonD[output_dir]}"

# 数据用于报告的链接目录，相对目录，不含子目录
link_dir <- "{config_jsonD[link_dir]}"

# 输出文件前缀
prefix <- "{config_jsonD[prefix]}"

# 根据前述定义自动生成的参数和文件名字，一般不需要修改
# 需要注意的是，原始count数据文件默认存储在数据目录的summary目录，后缀为.rc.xls

# 数据文件
cnt_file <- "{config_jsonD[cnt_file]}"
pheno_file <- "{config_jsonD[pheno_file]}"
feature_file <- "{config_jsonD[feature_file]}"

# 输出文件
rds_file <- paste0(output_dir, prefix, '.rds')
rds_filter_file <- paste0(output_dir, prefix, '.filter.rds')
cnt_filter_file <- paste0(output_dir,prefix,".filter_rc.xls")
pheno_filter_file <- paste0(output_dir,prefix,".filter_pheno.xls")

```

```{{r mkdir-output-dir}}
system(paste0("mkdir -p ", output_dir))
```

```{{r biology-common-data}}
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

样品信息`phenoData`通常包含以下信息 (Table \@ref(tab:pheno-example))。

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


```{{r load-data, cache.extra=mtime(c(cnt_file, pheno_file, feature_file))}}
if(debug) {{
  print(cnt_file)
}}
cnt_data <- read.table(cnt_file, sep="\\t", row.names=1, header=T,quote="", check.names=F)
cnt_data <- cnt_data[rowSums(cnt_data>0)>0,]

# Keep only samples appear in cnt_data and outut in same order as in cnt_data.
pheno_data <- read.table(pheno_file,sep="\\t", header=TRUE, row.names=1, quote="", check.names=F)
pheno_data <- pheno_data[match(colnames(cnt_data), rownames(pheno_data)),, drop=F]
pheno_data_col = colnames(pheno_data)
pheno_data_col_len <- length(pheno_data_col)
if(pheno_data_col_len<2) {{
  pheno_data_col <- rep(pheno_data_col,2)
}}

if((pheno_data_col[1]!="conditions") | (pheno_data_col[2]!="batch")) {{
  stop("***Pleace check if phenoData has conditions and batch at the second and third column.***")
}}


# Keep only features appeared in cnt_data and output in same order as in cnt_data.
feature_data <- read.table(feature_file, sep="\\t", header=T, row.names=1,quote="", check.names=F)
feature_data_rowname <- rownames(feature_data)
feature_data <- feature_data[match(rownames(cnt_data), feature_data_rowname),, drop=F]
feature_data_col <- colnames(feature_data)
if(! "Associated_Gene_Name" %in% feature_data_col){{
  stop("***Must have a column containing a list of unique gene symbols with column name as <Associated_Gene_Name>.***")
}}
```


基因表达Count数据如 (Table \@ref(tab:expr-cnt-part)) 所示。

```{{r expr-cnt-part}}
knitr::kable(
  head(cnt_data[,1:7]), booktabs=TRUE,
  caption="基因表达Count值展示 (前6行，前7列)。"
)
```


样品信息、处理方式、测序信息等如 (Table \@ref(tab:pheno-part))所示。

```{{r pheno-part}}
knitr::kable(
  head(pheno_data), booktabs=TRUE,
  caption="样品信息、处理方式、测序信息等 (前6行)。"
)
```

```{{asis echo=debug}}
把基因表达数据、样品信息数据、基因注释数据整合生成`SCESet`数据集。
```

```{{r sceSet}}
pheno_dataframe <- new("AnnotatedDataFrame", pheno_data)
feature_dataframe <- new("AnnotatedDataFrame", feature_data)

sceset_data <- scater::newSCESet(
  countData = cnt_data,
  phenoData = pheno_dataframe,
  featureData = feature_dataframe
)
```

移除在所有样品中表达值都为0的基因，获得用于下游分析的数据集如 (Table \@ref(tab:remove-all-zero))所示。

```{{r remove-all-zero}}
#keep_feature <- rowSums(counts(sceset_data)>0)>0
#sceset_data <- sceset_data[keep_feature,]
sta <- as.data.frame(dim(sceset_data))
rownames(sta) <- c("Number of genes","Number of samples")
colnames(sta) <- NULL
knitr::kable(t(sta), booktabs=TRUE,
  caption="表达矩阵统计。"
)
```

```{{asis echo=debug}}
定义对照基因集 - ERCC spikie-in 和线粒体基因。
```

```{{r ercc-mt}}
if (ercc) {{
  ERCC <- featureNames(sceset_data)[grepl("ERCC-", featureNames(sceset_data))]
  feature_controls = list(ERCC = ERCC, MT=mtgene)
}} else {{
  feature_controls = list(MT=mtgene)
}}
sceset_data <- scater::calculateQCMetrics(
    sceset_data,
    feature_controls = feature_controls)  
```

```{{asis echo=debug}}
查看sceset_data中包含的质量控制信息
```

```{{r sceset-data-qc-metric-list}}
feature_pheno = colnames(pData(sceset_data))
## [1] "total_counts"             "log10_total_counts"       "filter_on_total_counts"  
## [4] "total_features"           "log10_total_features"     "filter_on_total_features"
```

## 样品过滤 {{#sample-filter}}

样品在实验处理环节或者测序过程中存在一些不可预知的问题，从而产生一些异常的样品。为了保证后续结果的可靠性，需要通过判断样品测序文库大小、检测到的基因数目、ERCC spike-in基因的比例 (如果有的话)、线粒体基因的比例进行多层次的过滤。

通常使用中位数绝对偏差作为判断样品异常的标准。以测序文库大小为例，假设样品$i$中的Total read count是$r_i$，所有样品中Total read count的中位数是$m$，那么样品$i$ Total read count的绝对偏差$d_i$就是$d_{{i}}= |r_{{i}}-m|$。$d_{{i}}<`r nmads_seq_depth`*median(d_{{i}})$ 的样品会被移除 （移除测序深度低的样品）。为了增强过滤的鲁棒性，依据`样品测序的文库大小`和`检测到的基因数目`过滤时会先对相应对数值进行对数转换。依据`ERCC spike-in基因的比例`和`线粒体基因的比例`过滤时，$d_{{i}}>3*median(d_{{i}})$的样品会被移除 （移除检测到的内源基因少的样品）。

Figure \@ref(fig:hist-total-reads) 展示了每个样品检测到的reads数的分布，也就是测序文库的大小（如果使用了UMI，图形展示的就就是RNA分子的总数）。异常低测序深度的样品可能是因为捕获效率较低导致，通常会被过滤掉，过滤后结果如 (Table \@ref(tab:filter-total-counts-table))所示。



```{{r hist-total-reads, fig.cap="Histogram of library sizes for all samples. Red line represents meidan value of library sizes."}}
hist(sceset_data$total_counts/1e6, breaks=100, xlab="Library size of each sample (unit: million)", main="")
abline(v = median(sceset_data$total_counts/1e6), col = "red")
```

```{{r filter-total-counts-table}}
filter_by_library_size <- ! isOutlier(sceset_data$total_counts, nmads=nmads_seq_depth, type="lower", log=TRUE)
filter_by_library_size_table <- as.data.frame(table(filter_by_library_size))
filter_by_library_size_table$filter_by_library_size = revalue(filter_by_library_size_table$filter_by_library_size, c("FALSE"="Failed","TRUE"="Passed"), warn_missing = F)
colnames(filter_by_library_size_table) <- c("Library size filter", "Count of samples")
knitr::kable(
  filter_by_library_size_table,
  booktabs = TRUE,
  row.names = FALSE, format="html",
  caption = "The number of cells passed or failed the library size filter."
)
```

每个样品检测到的基因数目的分布如 (Figure \@ref(fig:hist-total-genes))所示。同文库大小信息类似，也可以用来过滤低质量样品, 结果如 (Table \@ref(tab:filter-by-expr-features-table))所示。

```{{r hist-total-genes, fig.cap="Histogram of the number of detected genes in all samples. Red line represents the median value of gene counts."}}
hist(sceset_data$total_features, breaks=100, xlab="Number of genes in each sample", main="")
abline(v = median(sceset_data$total_features), col = "red")
```

```{{r filter-by-expr-features-table}}
filter_by_expr_features <- ! isOutlier(sceset_data$total_features, nmads=nmads_gene_num, type="lower", log=TRUE)
filter_by_expr_features_table <- as.data.frame(table(filter_by_expr_features))
filter_by_expr_features_table$filter_by_expr_features = revalue(filter_by_expr_features_table$filter_by_expr_features, c("FALSE"="Failed","TRUE"="Passed"), warn_missing = F)
colnames(filter_by_expr_features_table) <- c("Detected genes filter", "Count of samples")
knitr::kable(
  filter_by_expr_features_table,
  booktabs = TRUE,
  row.names = FALSE,
  caption = "The number of cells passed or failed by total detected genes filter."
)
```

如果样品中有ERCC，就可以通过计算ERCC spike-in RNA和内源性RNA的比例评价细胞的质量。如果细胞中含有较高水平的ERCC spike-in RNA,那么这个细胞就会有较低的RNA起始量。这可能是因为细胞死亡或受到外界压力而导致了RNA降解，也有可能是捕获效率低。

一般来讲，如果样品中有高比例的reads比对到线粒体基因组，则意味着细胞出现了坏死或细胞质RNA丢失，这样的样品一般质量较差。但在使用线粒体基因占有的reads数过滤时需要充分考虑细胞的状态。如果细胞中有更多的线粒体基因或线粒体基因活性越强，其线粒体基因表达所占比例就会比较大。

假设大部分样品的建库质量都比较高并且样品状态也比较一致，有着过高或过低`ERCC spike-in RNA`和`线粒体RNA`的样品则被作为异常值去除，过滤标准分别是`r  nmads_ercc_percent`倍和`r nmads_mito_percent`倍的中值绝对偏差。

Figure \@ref(fig:pct-counts-feature-controls-MT-batch) 展示的是细胞中线粒体基因的表达量占细胞所有基因表达量的比例。Table \@ref(tab:filter-by-mt-table)展示的是根据异常线粒体基因表达量过滤掉的细胞。

分别使用不同的样品信息对样品进行颜色标记，用来检测异常细胞是否来源于同一批次或具有共同的样品特征。

```{{r pca-shapes-original}}
level <- length(unique(unlist(pData(sceset_data)[pheno_data_col[2]])))
shapes = (1:level)%%30  # maximum allow 30 types of symbols
# add following line to ggplot2 plot
#+ scale_shape_manual(values=shapes)
```

```{{r pct-counts-feature-controls-MT-batch, fig.height=9, fig.width=9, fig.cap="Percentage of read counts for mitochondrial genes relative to all genes. Each point represents a sample. The density plot at each axis represents the enrichment of samples. Samples are colored by different sample info to see which type of cells may be outliers."}}
#pheno_type <- "conditions"
scater::plotPhenoData(
    sceset_data,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_MT",
               colour = pheno_data_col[1],
               shape = pheno_data_col[2])
) + scale_shape_manual(values=shapes) + xlab("Number of detected genes") + ylab("Percentage of read counts for mitochondrial genes (%)") + guides(fill=guide_legend(ncol=2))
```



```{{r filter-by-mt-table}}
filter_by_MT <- ! isOutlier(sceset_data$pct_counts_feature_controls_MT, nmads=nmads_mito_percent, type="higher")
filter_by_MT_table <- as.data.frame(table(filter_by_MT))
filter_by_MT_table$filter_by_MT = revalue(filter_by_MT_table$filter_by_MT, c("FALSE"="Failed"
,"TRUE"="Passed"), warn_missing = F)
colnames(filter_by_MT_table) <- c("Mitochondrial gene filter", "Count of samples")
knitr::kable(
  filter_by_MT_table,
  booktabs = TRUE,
  row.names = FALSE,
  caption = "The number of cells passed or failed the mitochondrial gene filter. Cells with extremely high mitochonridal gene counts (as described above) will be removed."
)
```

```{{asis, echo=ercc}}
如果有ERCC spike-in, 可以评估覆盖ERCC spike-in 的read counts数占样品文库大小的比例。因为ERCC spike-in是等量加入的测序总RNA里面的，对于小的细胞或RNA降解严重的细胞中，其自身RNA含量少，ERCC read counts的比例会比较高。

Figure \@ref(fig:pct-counts-feature-controls-ERCC-FACS) 展示的是细胞中ERCC spike-in RNA
的表达量占细胞所有基因表达量的比例。

分别使用不同的样品信息对样品进行颜色标记，用来查看异常细胞是否来源于同一批次或具有共同的样品特征。
```

```{{r pct-counts-feature-controls-ERCC-FACS, fig.cap="Percentage of read counts for ERCC spike-in RNAs relative to all genes. Each point represents a sample. The density plot at each axis represents the enrichment of samples. Samples are colored by different sample info to see which type of cells may be outliers.", eval=ercc}}
#pheno_type <- "FACS"
scater::plotPhenoData(
    sceset_data,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_ERCC",
               colour = pheno_data_col[1],
               shape = pheno_data_col[2])
) + scale_shape_manual(values=shapes) + xlab("Number of detected genes") + ylab("Percentage of read counts for ERCC spike-in RNAs (%)")
```

```{{r filter-by-ercc-table, eval=ercc}}
filter_by_ERCC <- ! isOutlier(sceset_data$pct_counts_feature_controls_ERCC, nmads=nmads_ercc_percent, type="higher")
filter_by_ERCC_table <- as.data.frame(table(filter_by_ERCC))
filter_by_ERCC_table$filter_by_ERCC = revalue(filter_by_ERCC_table$filter_by_ERCC, c("FALSE"="Failed"
,"TRUE"="Passed"), warn_missing = F)
colnames(filter_by_ERCC_table) <- c("Mitochondrial gene filter", "Count of samples")
knitr::kable(
  filter_by_ERCC_table,
  booktabs = TRUE,
  row.names = FALSE,
  caption = "The number of cells passed or failed the mitochondrial gene filter. Cells with extremely high mitochonridal gene counts (as described above) will be removed."
)
```




```{{r define-filter}}
if (! ercc) {{
  filter_by_ERCC <- rep(TRUE, length(filter_by_expr_features))
}}
  sceset_data$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    
    # sufficient molecules counted (read counts)
    # keep samples with toral read counts larger than 80000
    filter_by_library_size &
  
    # remove cells with unusual number of reads in MT genes
    filter_by_MT &
    
    # sufficient endogenous RNA
    filter_by_ERCC
  )
```

Figure \@ref(fig:pca-outlier-filter) 展示了利用PCA自动识别出的异常样品与基于`样品测序的文库大小`、`检测到的基因数目`、`ERCC spike-in基因的比例`、`线粒体基因的比例`等识别出的异常样品的的吻合度。每个圆球状的点代表采用中位数绝对偏差过滤掉的样品，三角形点代表保留下的样品；橘黄色的点代表PCA自动检测出的异常点，蓝色的点代表通过PCA检测的样品。若PCA检测无异常样品，所有点颜色为黑色。通常情况下自动识别出的异常样品数目低于前面使用的过滤方法。


```{{asis, echo=FALSE}}

### Automatic threasholds

使用PCA自动识别潜在的异常细胞，默认情况下PCA分析会使用以细胞为行，下述变量为列的质量值矩阵，以2D图的形式展示。

* pct_counts_top_100_features
* total_features
* pct_counts_feature_controls
* n_detected_feature_controls  
* log10_counts_endogenous_features
* log10_counts_feature_controls
```

```{{r pca-outlier-filter, fig.cap="Sample outliers detected by MAD-based methods (use, FALSE (circles) represents filtered samples) and PCA-based methods (outlier, TRUE (orange) represents PCA deteced outliers)."}}
#selected_variables <- colnames(pData(sceset_data))[sapply(pData(sceset_data), is.numeric)]
#selected_variables <- selected_variables[selected_variables!="n_detected_feature_controls"]
selected_variables <- c("pct_counts_top_100_features","total_features","pct_counts_feature_controls","log10_counts_endogenous_features", "log10_counts_feature_controls")
sceset_data <- scater::plotPCA(sceset_data, 
                                          size_by="total_features",
                                          shape_by="use",
                                          pca_data_input="pdata",
                                          selected_variables=selected_variables,
                                          detect_outliers=TRUE,
                                          return_SCESet=TRUE)
filter_by_PCA <- sceset_data$outlier
```

```{{r define-filter2}}
sceset_data$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    
  # sufficient molecules counted (read counts)
  # keep samples with toral read counts larger than 80000
  filter_by_library_size &
  
  # remove cells with unusual number of reads in MT genes
  filter_by_MT &
    
  # sufficient endogenous RNA
  filter_by_ERCC &
    
  # PCA detected outliers
  (! filter_by_PCA)
)
```

Table \@ref(tab:filter-table) 展示不同的过滤标准分别过滤到的样品数目和最终用于后续分析的样品数目 （样品采用并集方式，未通过任意一个过滤标准的样品都不会保留）。

```{{r filter-table}}
#sceset_data <- sceset_data[, pData(sceset_data)$use]

filter_table = data.frame(Filtered_By_Lib_Size=sum(!filter_by_library_size), 
                                Filtered_By_Gene_Number=sum(!filter_by_expr_features),
                                Filtered_By_Mitochondrial_Gene_Percentage=sum(!filter_by_MT), 
                                Filtered_By_Spike_in_percentage=sum(!filter_by_ERCC), 
                                Filtered_By_PCA_autodetect=sum(filter_by_PCA),
                                Kept_Samples=ncol(sceset_data[, pData(sceset_data)$use]))
knitr::kable(
  t(filter_table),
  booktabs = TRUE,
  caption = 'Summary of filtered and kept samples.'
)
```


## 基因过滤

在移除低质量细胞的同时，也需要排除由于技术问题引入的异常表达的基因。
Figure \@ref(fig:top-50-gene-reads-count) 展示的是最高表达的50个基因的reads count的分布。如果单个基因表达量占总reads数比例低于1%-5%，则说明测序对全转录组的覆盖比较好。如果某一样品中最高表达基因的reads count占到total counts的20-50%，则说明这一样品文库质量较差导致测序存在偏好性。

此外，查看基因表达图谱也可以指导实验程序是否需要调整。
如果ERCC spike-in也出现在top 15表达的基因里面，说明加入的spike-in的量可以进一步稀释，既节约试剂，又能增大内源表达的RNA的检出量。


```{{r top-50-gene-reads-count, fig.cap="Number of total counts consumed by the top 50 expressed genes", fig.height=7}}
#如何显示基因名字，而不是ENSEMBLE编号
sceset_data_symbol <- sceset_data[, pData(sceset_data)$use]
featureNames(sceset_data_symbol) <- featureData(sceset_data_symbol)$Associated_Gene_Name
scater::plotQC(sceset_data_symbol, type = "highest-expression")
```

```{{r exprs-freq-vs-mean, fig.cap="Expression frequency showing number of cells with expression for the gene above zero or against mean expresssion level of the gene."}}
scater::plotQC(sceset_data, type="exprs-freq-vs-mean")
```


**通常还需要过滤掉表达水平极低的基因。**

如果使用了UMI，一般定义在至少两个细胞中检测到多于1个转录本并且平均转录本数目大于等于0.2的基因为表达的基因。

如果使用reads counts作为判断标准，定义在至少5个细胞中检测到多于1个reads并且平均read counts>1的基因为表达的基因。但这两种设定标准的方式都依赖于测序深度。

此外，要先过滤掉低质量的细胞，然后过滤不表达的基因，这样可以避免保留那些只在低质量的细胞中表达的基因。

过滤后保留和移除的基因如 Table \@ref(tab:gene-filter-reads)所示。

```{{r gene-filter-reads}}
if (UMI) {{
  read_count = 2
  mean_count = 0.2
}} else {{
  read_count = 5
  mean_count = 1
}}

filter_genes <- apply(counts(sceset_data[, pData(sceset_data)$use]),1,
                      function(x) ((length(x[x>1])>=read_count) & (mean(x)>=mean_count)))

fData(sceset_data)$use <- filter_genes

filter_genes_table <- as.data.frame(table(filter_genes))
filter_genes_table$filter_genes = revalue(filter_genes_table$filter_genes, c("FALSE"="Failed","TRUE"="Passed"), warn_missing = F)
colnames(filter_genes_table) <- c("Gene expression filter", "Count of genes")
knitr::kable(
  filter_genes_table,
  booktabs = TRUE,
  row.names = FALSE,
  caption = "The number of genes passed or failed gene filter."
)
```


```{{r gene-filter-umi, eval=UMI}}
#filter_genes <- apply(counts(sceset_data[, pData(sceset_data)$use]),1,
#                      function(x) ((length(x[x>1])>=2)

#fData(sceset_data)$use <- filter_genes

#knitr::kable(
#  as.data.frame(table(filter_genes)),
#  booktabs = TRUE,
#  row.names = FALSE,
#  caption = "'The number of genes removed by gene filter (FALSE)'"
#)
```

```{{asis, echo=FALSE}}
## 保存数据
```

```{{r save-data}}
saveRDS(sceset_data, file=rds_file)
```

```{{r filtered-data}}
a = dim(sceset_data[fData(sceset_data)$use, pData(sceset_data)$use])
geneNumber = as.vector(a[1]) 
sampleNumber = as.vector(a[2])
sceset_data_qc <- sceset_data[fData(sceset_data)$use, pData(sceset_data)$use]
saveRDS(sceset_data_qc, file=rds_filter_file)
sceset_data_qc_bak <- sceset_data[fData(sceset_data)$use, pData(sceset_data)$use]
if(ercc) {{
  endog_genes <- !fData(sceset_data_qc)$is_feature_control_ERCC
}} else {{
  endog_genes <- fData(sceset_data_qc)$use
}}

write.table(counts(sceset_data_qc[endog_genes,]), file=cnt_filter_file, sep="\\t", row.names = T, col.names=T, quote=F)
write.table(pData(sceset_data_qc)[,1:pheno_data_col_len],file=pheno_filter_file, sep="\\t",row.names=T,col.names=T,quote=F)
```

```{{r filtered-data-for-seurat}}
filtered_raw_count = counts(sceset_data_qc[endog_genes,])
filtered_raw_count_colname = colnames(filtered_raw_count)
conditions_filtered_raw_count = pData(sceset_data_qc[endog_genes, ])[[pheno_data_col[1]]]
filtered_raw_count_colname_seurat = paste(filtered_raw_count_colname, conditions_filtered_raw_count, sep="___")
filtered_raw_count_seurat <- filtered_raw_count
colnames(filtered_raw_count_seurat) <- filtered_raw_count_colname_seurat
rownames(filtered_raw_count_seurat) <- featureData(sceset_data_qc[endog_genes,])$Associated_Gene_Name
write.table(filtered_raw_count_seurat, file=paste0(output_dir,prefix,".seurat.filtered_rc.xls"),sep="\\t",row.names=T,col.names=T,quote=F)
seurat_pheno_filtered_raw_count = pData(sceset_data_qc[endog_genes, ])[,1:pheno_data_col_len]
seurat_pheno_filtered_raw_count$id = rownames(seurat_pheno_filtered_raw_count)
rownames(seurat_pheno_filtered_raw_count) <- filtered_raw_count_colname_seurat
write.table(seurat_pheno_filtered_raw_count, file=paste0(output_dir,prefix,".seurat.filtered_rc.pheno.xls"),sep="\\t",row.names=T,col.names=T,quote=F)

```

请点击下载过滤后基因的[Reads count](`r I(paste0(link_dir,basename(cnt_filter_file)))`)和样品的[pheno data](`r  I(paste0(link_dir,basename(pheno_filter_file)))`)信息。统计信息如 Table \@ref(tab:table-filtered-data-sum)。

```{{r table-filtered-data-sum}}
final_df <- data.frame("Kept samples"=sampleNumber, "Kept genes"=geneNumber)
rownames(final_df) <- "Count"
knitr::kable(final_df, booktabs=T, caption="Number of samples and genes passed all filters and used for downstream analysis.")
```


## 识别干扰因素

### 批次效应鉴定

单细胞测序数据分析的一个重要方面是控制数据的批次效应。批次效应是样品操作过程中引入的技术偏差，比如说，不同实验人员准备的两组样品，或不同的日期准备的两组样品，或不同的测序上机批次的样品，有时会发现同一时间或同一人员准备的样品具有更大的相似性，这会干扰真正的生物学变化。

因此要求在准备样品时详细记录收集样品、核酸提取、建库、测序等的时间，标记清楚哪些样品来源于同一批次，以便后期利用计算的方式屏蔽掉批次效应的影响，并且指导谨慎的解释结果。

```{{r read-in-data-qc}}
# 必要的时候读入数据，减少预处理时间
#sceset_data <- readRDS(paste0(dir, prefix,'.rds'))

# make clear use or use_default
#sceset_data_qc <- sceset_data[fData(sceset_data)$use, pData(sceset_data)$use]
#a = dim(sceset_data[fData(sceset_data)$use, pData(sceset_data)$use])
#geneNumber = as.vector(a[1]) 
#sampleNumber = as.vector(a[2])
#sceset_data_qc <- sceset_data[fData(sceset_data)$use, pData(sceset_data)$use]
#if(ercc) {{
#  endog_genes <- !fData(sceset_data_qc)$is_feature_control_ERCC
#}} else {{
#  endog_genes <- fData(sceset_data_qc)$use
#}}
```

可视化批次效应的一个方式是对样品进行聚类分析，其中一个方法是使用主成分分析 (PCA)。

PCA是一个统计分析方法，通过对原始数据的正交变换，把原始观察（基因表达）数据转换为一组线性无关的变量，称为主成分。每个主成分都能在一定程度上解释样品之间的差异。在做可视化时，一般选取第一、二两个主成分做二维散点图展示。具体参考 [GeneSino PCA](http://blog.genesino.com/2016/10/PCA/)。

Figure \@ref(fig:pca-after-qc-fig) 展示了完成质量控制用于下游分析的样品之间的相似性关系。每个点代表一个样品，点的大小代表样品中检测到的基因的数目；点的颜色和形状则代表样品不同的属性。(这个图基于基因的原始Reads count数生成。)

```{{r pca-shapes-qc}}
level <- length(unique(unlist(pData(sceset_data_qc)[pheno_data_col[2]])))
shapes = (1:level)%%30  # maximum allow 30 types of symbols
count_pheno_1 = length(unique(sceset_data_qc[[pheno_data_col[1]]]))
count_pheno_2 = length(unique(sceset_data_qc[[pheno_data_col[2]]]))
# add following line to ggplot2 plot
#+ scale_shape_manual(values=shapes)
```

```{{r pca-after-qc-fig, fig.height=11, fig.width=11, fig.cap="PCA plot of filtered data for showing batch effect."}}
if(count_pheno_1<10 | count_pheno_2<10){{
  if (count_pheno_2<10) {{
    colour_by=pheno_data_col[1]
    shape_by=pheno_data_col[2]
  }} else {{
    colour_by=pheno_data_col[2]
    shape_by=pheno_data_col[1]    
  }}
  scater::plotPCA(sceset_data_qc[endog_genes,], 
                ntop=5000, ncomponents = 3, 
                colour_by=colour_by,
                size_by="total_features",
                shape_by=shape_by,
                exprs_values="counts") + scale_shape_manual(values=shapes)
}} else {{
  scater::plotPCA(sceset_data_qc[endog_genes,], 
                ntop=5000, ncomponents = 3, 
                colour_by=pheno_data_col[1],
                size_by="total_features",
                exprs_values="counts")  
  scater::plotPCA(sceset_data_qc[endog_genes,], 
                ntop=5000, ncomponents = 3, 
                colour_by=pheno_data_col[2],
                size_by="total_features",
                exprs_values="counts")  
}} 

```




```{{asis echo=FALSE}}
### tSNE map 细胞分型

tSNE (t-Distributed Stochastic Neighbor Embedding) 通过组合降维（比如使用PCA）和最近邻网络随机行走策略在保留样品之间差异的情况下把高维数据（数万个基因的表达数据）映射到二维空间。tSNE是一个随机算法，每次运行结果都会有些差异，一般采用固定随机数发生器的种子来保证每次运行得到的结果一致。

`perplexity`用于设定邻近的点的数目，数值越小分出的类越多、网络越稀疏；数值越大网络越密集。最合适的`perplexity`值取决于数据的密度，数据密度越大，`perplexity`的值越大。具体参考 [Trickey TSNE](http://distill.pub/2016/misread-tsne/)。

Figure \@ref(fig:tsne-before-qc-fig) 展示了根据样品数目自动选取的多个`perplexity`获得的tSNE分类图。(这个图基于基因的原始Reads count数生成。)


Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
```

```{{r tsne-before-qc-fig, fig.cap="tSNE map of filtered data for different perplexity values.", eval=FALSE}}
#if (sampleNumber>20){{
#  count_perplexity = 20
#}} else {{
#  count_perplexity = sampleNumber
#}}
#cluster_perplexity = seq(5, round(count_perplexity/2),2)
#cluster_perplexity = round(sampleNumber / cluster_perplexity)
#cluster_perplexity = data.frame(a=cluster_perplexity)
#plot_tsne_list <- by(data=cluster_perplexity, INDICES = cluster_perplexity$a, FUN= function(x) {{
#  x <- x$a
#  a <- scater::plotTSNE(sceset_data_qc[endog_genes,], 
#           ntop=5000, perplexity = x,
#           colour_by=pheno_data_col[1],
#           shape_by=pheno_data_col[2],
#           exprs_values="counts",
#           rand_seed = 11521) + labs(title=paste0("perplexity = ",x))
#}})
#
#do.call(grid.arrange, c(plot_tsne_list, ncol=2))
```

### 干扰因素度量

单细胞测序过程中不可避免会引入一些干扰因素 (confounders)、人为操作误差和仪器偏好性等。由于每次操作的样品不可重复性，使得区分生物样品自身的差异和技术带来的偏差变得困难。


```{{asis, echo=FALSE}}
### Deteced genes

绘制检测的基因数与主成分的关系，可以看到基因数对PC1的贡献一致度在0.39。一般情况下，对于单细胞测序，检测到的基因数对于第一主成分的贡献都会很大，有的一致性可以达到0.8。
```


```{{r PC-correlation-detected-genes, eval=FALSE}}
scater::plotQC(sceset_data_qc[endog_genes,],
               type = "find-pcs",
               variable = "total_features",
               exprs_values = "counts")
```

```{{asis echo=FALSE}}
scater can also compute the marginal $R^2$ for each variable when fitting a linear model regressing expression values for each gene against just that variable, and display a density plot of the gene-wise marginal $R^2$ values for the variables.
```

因此需要对每个基因拟合一个线性回归模型反映不同的变量因素如检测到的基因总数、测序深度、测序批次对基因表达的贡献度来识别并移除这些干扰因素，指导数据的标准化和下游分析。

```{{asis echo=FALSE}}
This analysis indicates that the number of detected genes (again) and also the sequencing depth (number of counts) have substantial explanatory power for many genes, so these variables are good candidates for conditioning out in a normalisation step, or including in downstream statistical models. Expression of ERCCs also appears to be an important explanatory variable.
```

从Figure \@ref(fig:explanatory-variable) 可以看到，conditions (不同的样品)和Type（不同抗体富集）对基因差异的解释最大，这是符合预期的生物学本质差异；Total_counts、Total_features和Batch对基因差异的解释比较大，需要在标准化过程中考虑到。


```{{r explanatory-variable, fig.cap="Explanatory variables. Variables with peaks at the most right regions contribute the most to gene expression difference."}}
qc_variables <- c("total_features", "total_counts", unique(pheno_data_col), feature_pheno[grep("^pct_counts_feature_controls_", feature_pheno)])
#if (ercc) {{
#  qc_variables <- c(qc_variables, "pct_counts_technical_feature_controls")
#}}

scater::plotQC(sceset_data_qc[endog_genes,],
               type = "expl",
               exprs_values = "counts",
               variables = qc_variables) + ylab("Density")
```

```{{r plotExplanatoryVariables}}
#scater::plotExplanatoryVariables(sceset_data_qc)
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


