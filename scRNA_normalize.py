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
    This is designed to generate Rmd file for doing SC3 analysis.
'''

import sys
import os
from json import dumps as json_dumps
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
    parser.add_option("-r", "--rds-file", dest="rds_file",
        help="Filtered RDS file. Outlier samples and genes removed, ercc kept.")
    parser.add_option("-M", "--method", dest="method",
        help="scran or cpm")
    parser.add_option("-l", "--label", dest="label",
        default=' ', help="A unique string to represent the project. Only alphabets, numbers are allowed.")
    parser.add_option("-R", "--removeBatchEffect", dest="removeBatchEffect",
        type='int',
        help="whether you want to remove batch effect (1) or not (0).")
    #parser.add_option("-c", "--cluster-type", dest="cluster_type",
    #    help="<original> meaining show plots and find markers for original cluster. <recluster> meaning define new clusters.")
    #parser.add_option("-p", "--pheno-file", dest="pheno_file",
    #    help="(lowercase p) A pheno file according to format specified above. OPTIONAL.")
    parser.add_option("-o", "--output-dir", dest="output_directory",
        help="Specify the output directory of result file.")
    parser.add_option("-L", "--link-dir", dest="link_directory",
        help="Specify the linking directory of result file.")
    parser.add_option("-P", "--prefix", dest="prefix",
        help="(uppercase P) Prefix for output files.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.rds_file != None, "A label needed for -r"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    method = options.method
    removeBatchEffect = options.removeBatchEffect
    output_directory = options.output_directory
    link_directory   = options.link_directory
    prefix   = options.prefix
    rds_file = options.rds_file
    label    = options.label
    label_stub = label.replace(' ', '-')
    label_underline = label.replace(' ', '_')
    #matrix_file = options.matrix_file
    #pheno_file = options.pheno_file
    #log2_transform = options.log2_transform
    #cluster_type   = options.cluster_type
    #run_original = run_recluster = 'FALSE'
    #if cluster_type in ["original"]:
    #    run_original = 'TRUE'
    #elif cluster_type in ["recluster"]:
    #    run_recluster = 'TRUE'
    #else:
    #    print >>sys.stderr, "Unknown parameter for -c"

    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    print '''# {method} normalization

```{{r}}

if(! exists("output_dir")){{
    output_dir <- "{output_directory}"
}}

if(! exists("link_dir")){{
    link_dir <- "{link_directory}"
}}

if(! exists("prefix")){{
    prefix <- "{prefix}"
}}
```

```{{r initiate-{method}{label_stub}-data, cache.extra=mtime("{rds_file}")}}
sceset_data_{method}{label_underline}_qc <- readRDS("{rds_file}")

a_{method} = dim(sceset_data_{method}{label_underline}_qc)
geneNumber_{method} = as.vector(a_{method}[1]) 
sampleNumber_{method} = as.vector(a_{method}[2])

if(ercc) {{
  endog_genes_{method} <- !fData(sceset_data_{method}{label_underline}_qc)$is_feature_control_ERCC
}} else {{
  endog_genes_{method} <- fData(sceset_data_{method}{label_underline}_qc)$use
}}

## Uncomment when necessary
# level <- length(unique(unlist(pData(sceset_data_{method}{label_underline}_qc)[pheno_data_col[2]])))
# shapes = (1:level)%%30  # maximum allow 30 types of symbols
# add following line to ggplot2 plot
#+ scale_shape_manual(values=shapes)
```

```{{r table-filtered-data-sum-{method} }}
final_df_{method} <- data.frame("Kept samples"=sampleNumber_{method}, "Kept genes"=geneNumber_{method})
rownames(final_df_{method}) <- "Count"
knitr::kable(final_df_{method}, booktabs=T, caption="Number of samples and genes passed all filters and used for downstream analysis.")
```


```{{r}}
calc_cell_RLE <- function (expr_mat, spikes=NULL){{
  RLE_gene <- function(x) {{
    if (median(unlist(x)) > 0){{
      log((x+1)/(median(unlist(x))+1))/log(2)
    }} else {{
      rep(NA, times=length(x))
    }}
  }}
  if (!is.null(spikes)) {{
    # Each row (expression of one gene in different samples) divide by the median of this row then get log2 value
    # If genes are normalized into median expression, the final log2 value would be claso to zero.
    RLE_matrix <- t(apply(expr_mat[-spikes,],1,RLE_gene))
  }} else {{
    RLE_matrix <- t(apply(expr_mat,1,RLE_gene))
  }}
  # Get median value for each sample.
  # Theoritically the median for each sample should be all zero.
  cell_RLE <- apply(RLE_matrix, 2, median, na.rm=T)
  return(cell_RLE)
}}

calc_cell_RLE2 <- function (expr_mat, spikes=NULL){{
  RLE_gene <- function(x) {{
    if (median(unlist(x)) > 0){{
      log((x+1)/(median(unlist(x))+1))/log(2)
    }} else {{
      rep(NA, times=length(x))
    }}
  }}
  if (!is.null(spikes)) {{
    # Each row (expression of one gene in different samples) divide by the median of this row then get log2 value
    # If genes are normalized into median expression, the final log2 value would be claso to zero.
    RLE_matrix <- t(apply(expr_mat[-spikes,],1,RLE_gene))
  }} else {{
    RLE_matrix <- t(apply(expr_mat,1,RLE_gene))
  }}
  return(RLE_matrix)
}}
```
'''.format(method=method, rds_file=rds_file, label_stub=label_stub,  label_underline=label_underline, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix)

    if method == 'scran':
        print '''
# 基因表达标准化 ({method})

**{method}**是CPM (counts per million) 应用于单细胞测序分析的变体。这个方法通过把多个单细胞样品聚合到一起计算其总reads数作为归一化因子。因为一个细胞可以出现在多个pool里面，细胞类型特异的因子可以通过线性代数去卷积的方法从pool特异性的因子计算而来 [@L.Lun2016]。

我们使用**{method}**做文库标准化,并使用PCA plot (Figure \@ref(fig:pca-after-qc-fig-for-{method}))和计算标准化后的相对log表达 (Relative log expression RLE) (Figure \@ref(fig:cell-wise-RLE-{method}))来比较文库归一化的效果。如果样品测序了更多的reads，就会有更多的基因表达值高于中值，使得RLE值为正值；反之，如果样品测序了较低的深度，就会有更多的基因表达值低于中值，使得RLE为负；归一化的细胞的RLE值应接近0。


```{{asis echo=FALSE}}
We will use visual inspection of PCA plots and calculation of cell-wise relative log expression (calc_cell_RLE()) to compare the efficiency of different normalization methods. Cells with many[few] reads have higher[lower] than median expression for most genes resulting in a positive[negative] RLE across the cell, whereas normalized cells have an RLE close to zero.
```

```{{r qc-{method}{label_stub}-normalize}}
if (sampleNumber_{method} > 600){{
  min_size_each_cluster <- 60
  sizes <- c(10,20,30,40,50)
  qclust <- {method}::quickCluster(sceset_data_{method}{label_underline}_qc, min.size = min_size_each_cluster)
  sceset_data_{method}{label_underline}_qc <- {method}::computeSumFactors(sceset_data_{method}{label_underline}_qc, sizes=sizes, clusters=qclust,positive=TRUE)  
}} else if (sampleNumber_{method} > 100) {{
  min_size_each_cluster <- 40
  sizes <- c(10,20,30)
  qclust <- {method}::quickCluster(sceset_data_{method}{label_underline}_qc, min.size = min_size_each_cluster)
  sceset_data_{method}{label_underline}_qc <- {method}::computeSumFactors(sceset_data_{method}{label_underline}_qc, sizes=sizes, clusters=qclust,positive=TRUE)
}} else {{
  sizes = seq(5, sampleNumber_{method}, 10)
  sceset_data_{method}{label_underline}_qc <- {method}::computeSumFactors(sceset_data_{method}{label_underline}_qc, sizes=sizes,positive=TRUE)
}}

summary(sizeFactors(sceset_data_{method}{label_underline}_qc))

sizeFactorsScran <- as.data.frame(sizeFactors(sceset_data_{method}{label_underline}_qc))
colnames(sizeFactorsScran) <- c("size_factor_{method}")
write.table(sizeFactorsScran, file=paste0(output_dir,prefix,".size_factor_{method}.xls"), row.names=T,col.names=T,sep="\t",quote=F)


#sceset_data_{method}{label_underline}_qc <- {method}::computeSumFactors(sceset_data_{method}{label_underline}_qc, sizes=pool_size, clusters=qclust, positive=TRUE)

#sceset_data_{method}{label_underline}_qc <- {method}::computeSumFactors(sceset_data_{method}{label_underline}_qc, sizes=c(20,40,60,80), positive=TRUE)
if (ercc) {{
  sceset_data_{method}{label_underline}_qc <- computeSpikeFactors(sceset_data_{method}{label_underline}_qc, type="ERCC", general.use=FALSE)
}}

# The count data are used to compute normalized log-expression values for use in downstream analyses. Each value is defined as the log-ratio of each count to the size factor for the corresponding cell, after adding a prior count of 1 to avoid undefined values at zero counts. Division by the size factor ensures that any cell-specific biases are removed. If spike-in-specific size factors are present in sce, they will be automatically applied to normalize the spike-in transcripts separately from the endogenous genes.
#The log-transformation provides some measure of variance stabilization (Law et al. 2014), so that high-abundance genes with large variances do not dominate downstream analyses. The computed values are stored as an exprs matrix in addition to the other assay elements.
sceset_data_{method}{label_underline}_qc <- scater::normalize(sceset_data_{method}{label_underline}_qc)


#Error in svd(x, nu = 0) : a dimension is zero
```
'''.format(method=method, label_stub=label_stub,  label_underline=label_underline)
    elif method == "cpm":
        print '''

我们使用**CPM**做文库标准化,并使用PCA plot (Figure \@ref(fig:pca-after-qc-fig-for-cpm))和计算标准化后的相对log表达(Relative log expression RLE) (Figure \@ref(fig:cell-wise-RLE-cpm))来比较文库归一化的效果。如果样品测序了更多的reads，就会有更多的基因表达值高于中值，使得RLE值为正值；反之，如果样品测序了较低的深度，就会有更多的基因表达值低于中值，使得RLE为负；归一化的细胞的RLE值应接近0。

```{{asis echo=FALSE}}
We will use visual inspection of PCA plots and calculation of cell-wise relative log expression (calc_cell_RLE()) to compare the efficiency of different normalization methods. Cells with many[few] reads have higher[lower] than median expression for most genes resulting in a positive[negative] RLE across the cell, whereas normalized cells have an RLE close to zero.
```

```{{r ${method}{label_stub}-norm-log2}}
# Default logged, no this operation needed.
# set_exprs(sceset_data_{method}{label_underline}_qc, "exprs") = log2(exprs(sceset_data_{method}{label_underline}_qc)+1)
#sceset_data_{method}{label_underline}_qc@logged=TRUE
```

'''.format(method=method, label_stub=label_stub,  label_underline=label_underline)
    else:
        print >>sys.stderr, "Unknown parameters"
        sys.exit(1)

    print '''
```{{r save-rds-{method}{label_stub}-norm-endog}}
sceset_data_{method}{label_underline}_qc_endog <- sceset_data_{method}{label_underline}_qc[endog_genes_{method}, ]    
saveRDS(sceset_data_{method}{label_underline}_qc_endog, file=paste0(output_dir, prefix, ".{method}.norm.log2.endog.rds"))
```

```{{r pca-after-qc-fig-for-{method}, fig.height=11, fig.width=11, fig.cap="PCA plot of {method} normalized data."}}
if(count_pheno_1<10 | count_pheno_2<10){{
  if (count_pheno_2<10) {{
    colour_by=pheno_data_col[1]
    shape_by=pheno_data_col[2]
  }} else {{
    colour_by=pheno_data_col[2]
    shape_by=pheno_data_col[1]    
  }}
  scater::plotPCA(sceset_data_{method}{label_underline}_qc_endog, 
                ntop=5000, ncomponents = 4, 
                colour_by=colour_by,
                size_by="total_features",
                shape_by=shape_by,
                exprs_values="exprs") + scale_shape_manual(values=shapes) + guides(fill=guide_legend(ncol=2))
}} else {{
  scater::plotPCA(sceset_data_{method}{label_underline}_qc_endog, 
                ntop=5000, ncomponents = 4, 
                colour_by=pheno_data_col[1],
                size_by="total_features",
                exprs_values="exprs")  + guides(fill=guide_legend(ncol=2))
  scater::plotPCA(sceset_data_{method}{label_underline}_qc_endog, 
                ntop=5000, ncomponents = 4, 
                colour_by=pheno_data_col[2],
                size_by="total_features",
                exprs_values="exprs") + guides(fill=guide_legend(ncol=2)) 
}} 

```

理想条件下，标准化后的RLE的值应该为0。如果非0，说明标准化方法不合适或部分样品还存在干扰因素如Batch effect没有移除。

```{{r cell-wise-RLE-{method}, fig.cap="Cell-wise RLE of {method} normalized data. More zero, more bettter."}}
{method}_normalized_log2_expr = exprs(sceset_data_{method}{label_underline}_qc[endog_genes_{method}, ])

raw_rle = calc_cell_RLE(counts(sceset_data_{method}{label_underline}_qc[endog_genes_{method}, ]))
{method}_rle = calc_cell_RLE({method}_normalized_log2_expr)

rle_mat_{method} <- data.frame(variable=c(rep("Raw", length(raw_rle)), c(rep("Normalized", length({method}_rle)))), value=c(raw_rle, {method}_rle))
rle_mat_{method}$variable <- factor(rle_mat_{method}$variable, levels=c("Raw","Normalized"))
ggplot(rle_mat_{method}, aes(variable,value)) + ylab("Relative log expression (RLE).") + xlab("") + 
  theme_bw() + theme(panel.grid.major = element_blank()) +
  theme(legend.position="none") + geom_boxplot()
```

标准化并且log2转换的数据下载 [normalized_log2_transfered](`r I(paste0(link_dir,prefix,".{method}_normalized_log2transformed.xls"))`)。

```{{r output-normalized-data-{method}}}
write.table({method}_normalized_log2_expr, file=paste0(output_dir,prefix,".{method}_normalized_log2transformed.xls"), row.names=T,col.names=T,sep="\\t",quote=F)

{method}_normalized_log2_expr_colname = colnames({method}_normalized_log2_expr)
conditions_tmp_{method}{label_underline}_norm = pData(sceset_data_{method}{label_underline}_qc_endog)[[pheno_data_col[1]]]
{method}_normalized_log2_expr_colname_seurat = paste({method}_normalized_log2_expr_colname,conditions_tmp_{method}{label_underline}_norm,sep="___")
{method}_normalized_log2_expr_seurat <- {method}_normalized_log2_expr
colnames({method}_normalized_log2_expr_seurat) <- {method}_normalized_log2_expr_colname_seurat
rownames({method}_normalized_log2_expr_seurat) <- featureData(sceset_data_{method}{label_underline}_qc_endog)$Associated_Gene_Name
write.table({method}_normalized_log2_expr_seurat, file=paste0(output_dir,prefix,".seurat.{method}_norm_log2.xls"),sep="\\t",row.names=T,col.names=T,quote=F)
seurat_pheno_{method}{label_underline}_norm = pData(sceset_data_{method}{label_underline}_qc_endog)[,1:pheno_data_col_len]
seurat_pheno_{method}{label_underline}_norm$id = rownames(seurat_pheno_{method}{label_underline}_norm)
rownames(seurat_pheno_{method}{label_underline}_norm) <- {method}_normalized_log2_expr_colname_seurat
write.table(seurat_pheno_{method}{label_underline}_norm, file=paste0(output_dir,prefix,".seurat.{method}_norm_log2.pheno.xls"),sep="\\t",row.names=T,col.names=T,quote=F)
```

## 基因表达标准化后细胞预分型

tSNE (t-Distributed Stochastic Neighbor Embedding) 通过组合降维（比如使用PCA）和最近邻网络随机行走策略在保留样品之间差异的情况下把高维数据（数万个基因的表达数据）映射到二维空间。tSNE是一个随机算法，每次运行结果都会有些差异，一般采用固定随机数发生器的种子来保证每次运行得到的结果一致。

`perplexity`用于设定邻近的点的数目，数值越小分出的类越多、网络越稀疏；数值越大网络越密集。最合适的`perplexity`值取决于数据的密度，数据密度越大，`perplexity`的值越大。具体参考 [Trickey TSNE](http://distill.pub/2016/misread-tsne/)。

Figure \@ref(fig:tsne-cell-typing-normalized-{method}) 展示了根据样品数目自动选取的多个`perplexity`获得的tSNE分类图。(这个图基于{method} normalized的基因表达生成。)

```{{asis echo=FALSE}}
Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
```

```{{r cluster-perplexity-{method}}}
step_perplexity_{method} = round(sampleNumber_{method}/10)
if (step_perplexity_{method}<2){{
  step_perplexity_{method} = 2
}}
count_perplexity_{method} = round(sampleNumber_{method}/2)
count_midpoint_{method} = count_perplexity_{method}
if (count_midpoint_{method}>15){{
  count_midpoint_{method}=15
}}

cluster_perplexity_{method} = c(seq(5,count_midpoint_{method},3),seq(count_midpoint_{method}, count_perplexity_{method},step_perplexity_{method}))
cluster_perplexity_{method} = unique(round(sampleNumber_{method} / cluster_perplexity_{method}))
cluster_perplexity_{method} = data.frame(a=cluster_perplexity_{method})
```

```{{r tsne-cell-typing-normalized-{method}, fig.height=16, fig.cap="tSNE map of normalized expression data given different perplexity values."}}
plot_tsne_list_{method} <- by(data=cluster_perplexity_{method}, INDICES = cluster_perplexity_{method}$a, FUN= function(x) {{
  x <- x$a
  if(count_pheno_1<10 | count_pheno_2<10){{
    if (count_pheno_2<10) {{
      colour_by=pheno_data_col[1]
      shape_by=pheno_data_col[2]
    }} else {{
      colour_by=pheno_data_col[2]
      shape_by=pheno_data_col[1]    
    }}
    a <- scater::plotTSNE(sceset_data_{method}{label_underline}_qc_endog, 
         ntop=500, perplexity = x,
         colour_by=colour_by,
         shape_by=shape_by,
         exprs_values="exprs",
         rand_seed = 11521) + scale_shape_manual(values=shapes) + labs(title=paste0("perplexity = ",x)) + guides(fill=guide_legend(ncol=2))
  }} else {{
    a <- scater::plotTSNE(sceset_data_{method}{label_underline}_qc_endog, 
         ntop=500, perplexity = x,
         colour_by=pheno_data_col[1],
         exprs_values="exprs",
         rand_seed = 11521) + labs(title=paste0("perplexity = ",x)) + guides(fill=guide_legend(ncol=2))    
  }}
}})

do.call(grid.arrange, c(plot_tsne_list_{method}, ncol=2))

```

## 移除批次效应后细胞预分型

```{{r remove-batch-effect-{method}}}
batch_{method} = sceset_data_{method}{label_underline}_qc_endog[[pheno_data_col[2]]]
exprs_batch_corrected_{method} <- limma::removeBatchEffect({method}_normalized_log2_expr, batch=batch_{method})
#set_exprs(sceset_data_{method}{label_underline}_qc,"{method}_limma_batch_removed_log2") = exprs_batch_corrected_{method}
sceset_data_{method}{label_underline}_qc_batchR_endog <- sceset_data_{method}{label_underline}_qc_endog
set_exprs(sceset_data_{method}{label_underline}_qc_batchR_endog, "exprs") = exprs_batch_corrected_{method}

saveRDS(sceset_data_{method}{label_underline}_qc_batchR_endog, file=paste0(output_dir, prefix, ".{method}.norm.log2.batchRm.endog.rds"))

exprs_batch_corrected_{method}{label_underline}_colname = colnames(exprs_batch_corrected_{method})
conditions_tmp_{method}{label_underline}_norm_batchlmR = pData(sceset_data_{method}{label_underline}_qc_batchR_endog)[[pheno_data_col[1]]]
exprs_batch_corrected_{method}{label_underline}_colname_seurat = paste(exprs_batch_corrected_{method}{label_underline}_colname,conditions_tmp_{method}{label_underline}_norm_batchlmR,sep="___")
exprs_batch_corrected_{method}{label_underline}_seurat <- exprs_batch_corrected_{method}
colnames(exprs_batch_corrected_{method}{label_underline}_seurat) <- exprs_batch_corrected_{method}{label_underline}_colname_seurat
rownames(exprs_batch_corrected_{method}{label_underline}_seurat) <- featureData(sceset_data_{method}{label_underline}_qc_batchR_endog)$Associated_Gene_Name
write.table(exprs_batch_corrected_{method}{label_underline}_seurat, file=paste0(output_dir,prefix,".seurat.{method}_norm_log2_batchRm.xls"),sep="\\t",row.names=T,col.names=T,quote=F)
seurat_pheno_{method}{label_underline}_norm_batchlmR = pData(sceset_data_{method}{label_underline}_qc_batchR_endog)[,1:pheno_data_col_len]
seurat_pheno_{method}{label_underline}_norm_batchlmR$id = rownames(seurat_pheno_{method}{label_underline}_norm_batchlmR)
rownames(seurat_pheno_{method}{label_underline}_norm_batchlmR) <- exprs_batch_corrected_{method}{label_underline}_colname_seurat
write.table(seurat_pheno_{method}{label_underline}_norm_batchlmR, file=paste0(output_dir,prefix,".seurat.{method}_norm_log2_batchRm.pheno.xls"),sep="\\t",row.names=T,col.names=T,quote=F)
```

采用拟合分析方法，根据实验记录的样品批次，移除批次效应对数据的影响，获得标准化并且log转换并且移除批次效应的数据，可点击链接下载 [normalized_batch_removed_log2transfered](`r I(paste0(link_dir,prefix,".{method}_norm_log2_batchRm.xls"))`)。

```{{r exprs-batch-corrected-{method}}}
write.table(exprs_batch_corrected_{method}, file=paste0(output_dir,prefix,".{method}_norm_log2_batchRm.xls"),sep="\\t",row.names=T,col.names=T,quote=F)
```

Figure \@ref(fig:tsne-cell-typing-{method}{label_stub}-normalized-batch-removed) 展示了根据样品数目自动选取的多个`perplexity`获得的tSNE分类图。(这个图基于移除批次效应的标准化基因表达数据生成。)


```{{r tsne-cell-typing-{method}{label_stub}-normalized-batch-removed, fig.height=16, fig.cap="tSNE map of normalized expression data with batch effect removed given different perplexity values."}}

plot_tsne_list_{method} <- by(data=cluster_perplexity_{method}, INDICES = cluster_perplexity_{method}$a, FUN= function(x) {{
  x <- x$a
  if(count_pheno_1<10 | count_pheno_2<10){{
    if (count_pheno_2<10) {{
      colour_by=pheno_data_col[1]
      shape_by=pheno_data_col[2]
    }} else {{
      colour_by=pheno_data_col[2]
      shape_by=pheno_data_col[1]    
    }}
    a <- scater::plotTSNE(sceset_data_{method}{label_underline}_qc_batchR_endog, 
         ntop=500, perplexity = x,
         colour_by=colour_by,
         shape_by=shape_by,
         exprs_values="exprs",
         rand_seed = 11521) + scale_shape_manual(values=shapes) + labs(title=paste0("perplexity = ",x)) + guides(fill=guide_legend(ncol=2))
  }} else {{
    a <- scater::plotTSNE(sceset_data_{method}{label_underline}_qc_batchR_endog, 
         ntop=500, perplexity = x,
         colour_by=pheno_data_col[1],
         exprs_values="exprs",
         rand_seed = 11521) + labs(title=paste0("perplexity = ",x))+guides(fill=guide_legend(ncol=2))    
  }}
}})

do.call(grid.arrange, c(plot_tsne_list_{method}, ncol=2))

```{{r}}
system(paste0("mkdir -p ", link_dir))
system(paste0("/bin/cp -up ", output_dir,"* ", link_dir))
```

'''.format(rds_file=rds_file, method=method, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, label_stub=label_stub, label_underline=label_underline)
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


