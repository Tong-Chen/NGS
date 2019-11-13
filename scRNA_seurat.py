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
    This is designed to generate Rmd file for doing seurat analysis.
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
    parser.add_option("-l", "--label", dest="label",
        help="several concise words to represent the content of this analysis. Only alphabets and blanks are allowed. Must be unique if this script will be run multiple times. String like <scran norm original cluster>, <cpm norm reclassify>, <scran norm batch remove reclassify>,  <cpm norm sc3 cluster>")
    parser.add_option("-d", "--description", dest="descrip",
        help="A paragraph of markdown formated words to describe things one wants to do.")
    parser.add_option("-m", "--matrix-file", dest="matrix_file",
        help="A matrix file according to format specified above.")
    parser.add_option("-s", "--matrix-name-separtor", dest="matrix_delim",
        default='___', 
        help="The delim used to separate samples and their group information. \
This name attribute line locates in the first row of matrix file. Default <___>.")
    parser.add_option("-M", "--mito-genes", dest="mito_gene",
        help="A list of mitochondrial genes (gene name should match thos in matrix file, separated by comma.). or a string <^MT-> (the program will grep genes which name starts with MT-.)")
    parser.add_option("-e", "--max-allowed-genes", dest="maximum_gene",
        default=15000, type='int', help="Maximum allowed genes. Default 15000.")
    parser.add_option("-E", "--max-allowed-mito-percent", dest="max_mito_percent",
        default=0.7, type='float', help="Maximum allowed mito gene percent. Default 0.7.")
    parser.add_option("-p", "--pheno-file", dest="pheno_file",
        help="(lowercase p) A pheno file according to format specified above. OPTIONAL.")
    parser.add_option("-c", "--cluster-type", dest="cluster_type",
        help="<recluster> or <original>.")
    parser.add_option("-o", "--output-dir", dest="output_directory",
        help="Specify the output directory of result file.")
    parser.add_option("-L", "--link-dir", dest="link_directory",
        help="Specify the linking directory of result file.")
    parser.add_option("-P", "--prefix", dest="prefix",
        help="(uppercase P) Prefix for output files.")
    parser.add_option("-t", "--log2-transform", dest="log2_transform",
        help="Do the data need log2_transform. Accept TRUE to transform. Specfify FALSE if already transformed.")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-D", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.label != None, "A label needed for -l"
    assert options.log2_transform != None, "Parameter needed for -t"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    cluster_type   = options.cluster_type
    label = options.label + ' ' + cluster_type
    label_underline = label.replace(' ', '_')
    label_stub = label.replace(' ', '-') 
    descrip = options.descrip
    output_directory = options.output_directory
    link_directory   = options.link_directory
    prefix   = options.prefix
    matrix_file = options.matrix_file
    delim = options.matrix_delim
    mito_gene = options.mito_gene
    pheno_file = options.pheno_file
    log2_transform = options.log2_transform
    maximum_gene = options.maximum_gene 
    max_mito_percent = options.max_mito_percent  
    run_original = run_recluster = 'FALSE'
    if cluster_type in ["original"]:
        run_original = 'TRUE'
    elif cluster_type in ["recluster"]:
        run_recluster = 'TRUE'
    else:
        print >>sys.stderr, "Unknown parameter for -c"

    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    print '''
```{{r}}
library(Seurat)
library(Matrix)
library(dplyr)

if(! exists("output_dir")){{
    output_dir <- "{output_directory}"
    system("mkdir -p {output_directory}")
}}

if(! exists("link_dir")){{
    link_dir <- "{link_directory}"
    system("mkdir -p {link_directory}")
}}

if(! exists("prefix")){{
    prefix <- "{prefix}"
}}


```

# {descrip} ({label})

{descrip}

'''.format(label=label, label_underline=label_underline, label_stub=label_stub, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, matrix_file=matrix_file, pheno_file=pheno_file, 
        log2_transform=log2_transform, run_recluster=run_recluster, 
        run_original=run_original, descrip=descrip, cluster_type=cluster_type, delim=delim)
    
    if pheno_file:
        print ''' 

```{{r {label_stub}-data-seaurt,cache.extra=mtime(c("{matrix_file}","{pheno_file}"))}}
seurat_raw_{label_underline} <- read.table("{matrix_file}", sep="\\t", header=T, row.names=1, quote="")

seurat_pheno_{label_underline} <- read.table("{pheno_file}", sep="\\t", header=T, row.names=1, quote="")
seurat_samplenumber_{label_underline} <- dim(seurat_raw_{label_underline})[2]

step_perplexity_{label_underline} = round(seurat_samplenumber_{label_underline}/10)
if (step_perplexity_{label_underline}<2){{
  step_perplexity_{label_underline} = 2
}}
count_perplexity_{label_underline} = round(seurat_samplenumber_{label_underline}/2)
count_midpoint_{label_underline} = count_perplexity_{label_underline}
if (count_midpoint_{label_underline}>15){{
    count_midpoint_{label_underline}=15
}}
    
cluster_perplexity_{label_underline} = c(seq(5, count_midpoint_{label_underline}, 3), seq(count_midpoint_{label_underline},  count_perplexity_{label_underline}, step_perplexity_{label_underline}))
cluster_perplexity_{label_underline} = unique(round(seurat_samplenumber_{label_underline} / cluster_perplexity_{label_underline}))
cluster_perplexity_{label_underline} = data.frame(a=cluster_perplexity_{label_underline})

```


```{{r {label_stub}-seurat-initial, include=FALSE}}
seurat_ct_{label_underline} <- new('seurat', raw.data=seurat_raw_{label_underline})
seurat_ct_{label_underline} <- Seurat::Setup(seurat_ct_{label_underline}, project=prefix, min.cells=1, min.genes=1, names.field=2, names.delim="{delim}", is.expr=0, do.logNormalize={log2_transform}, meta.data=seurat_pheno_{label_underline}) 
#seurat_ct_{label_underline}_largestData <- max(max(seurat_ct_{label_underline}@scale.data))

```

'''.format(label=label, label_underline=label_underline, label_stub=label_stub, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, matrix_file=matrix_file, pheno_file=pheno_file, 
        log2_transform=log2_transform, run_recluster=run_recluster, 
        run_original=run_original, descrip=descrip, cluster_type=cluster_type, delim=delim)
    else:
        print ''' 

```{{r {label_stub}-data-seaurt,cache.extra=mtime(c("{matrix_file}"))}}
seurat_raw_{label_underline} <- read.table("{matrix_file}", sep="\\t", header=T, row.names=1, quote="")

seurat_samplenumber_{label_underline} <- dim(seurat_raw_{label_underline})[2]

step_perplexity_{label_underline} = round(seurat_samplenumber_{label_underline}/10)
if (step_perplexity_{label_underline}<2){{
  step_perplexity_{label_underline} = 2
}}
count_perplexity_{label_underline} = round(seurat_samplenumber_{label_underline}/2)
count_midpoint_{label_underline} = count_perplexity_{label_underline}
if (count_midpoint_{label_underline}>15){{
    count_midpoint_{label_underline}=15
}}
    
cluster_perplexity_{label_underline} = c(seq(5, count_midpoint_{label_underline}, 3), seq(count_midpoint_{label_underline},  count_perplexity_{label_underline}, step_perplexity_{label_underline}))
cluster_perplexity_{label_underline} = unique(round(seurat_samplenumber_{label_underline} / cluster_perplexity_{label_underline}))
cluster_perplexity_{label_underline} = data.frame(a=cluster_perplexity_{label_underline})

```


```{{r {label_stub}-seurat-initial, include=FALSE}}
seurat_ct_{label_underline} <- new('seurat', raw.data=seurat_raw_{label_underline})
seurat_ct_{label_underline} <- Seurat::Setup(seurat_ct_{label_underline}, project=prefix, min.cells=1, min.genes=1, names.field=2, names.delim="{delim}", is.expr=0, do.logNormalize={log2_transform}) 
#seurat_ct_{label_underline}_largestData <- max(max(seurat_ct_{label_underline}@scale.data))

```

'''.format(label=label, label_underline=label_underline, label_stub=label_stub, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, matrix_file=matrix_file, pheno_file=pheno_file, 
        log2_transform=log2_transform, run_recluster=run_recluster, 
        run_original=run_original, descrip=descrip, cluster_type=cluster_type, delim=delim)


    if mito_gene:
        if mito_gene == "^MT-":
            mito_gene = 'mito.genes <- grep("^MT-", rownames(seurat_ct_{label_underline}@data), value=T)'.format(label_underline=label_underline)
        else:
            mito_geneL = ','.join(["\""+i.strip()+"\"" for i in mito_gene.strip().split(',')])
            mito_gene = 'mito.genes <- c('+ mito_geneL +')'

        print '''

```{{r {label_stub}-mito-gene}}        
{mito_gene}
percent.mito <- colSums(expm1(seurat_ct_{label_underline}@data[mito.genes,]))/colSums(expm1(seurat_ct_{label_underline}@data))
seurat_ct_{label_underline} <- Seurat::AddMetaData(seurat_ct_{label_underline}, percent.mito, "percent.mito")
seurat_ct_{label_underline} <- Seurat::SubsetData(seurat_ct_{label_underline}, subset.name="nGene", accept.high={maximum_gene})
seurat_ct_{label_underline} <- Seurat::SubsetData(seurat_ct_{label_underline}, subset.name="percent.mito", accept.high={max_mito_percent})
#seurat_ct_{label_underline}2 <- Seurat::RegressOut(seurat_ct_{label_underline}, latent.vars=c("percent.mito"), do.scale=F, do.center=F)
#seurat_ct_{label_underline}2_regress <- as(seurat_ct_{label_underline}2@scale.data, matrix)
#write.table(seurat_ct_{label_underline}2_regress, file=paste0(output_dir,prefix,".{label_underline}.seurat_{cluster_type}_mito_regress_no_scale_expr.xls"), sep="\\t", row.names=T, col.names=T, quote=F)
seurat_ct_{label_underline} <- Seurat::RegressOut(seurat_ct_{label_underline}, latent.vars=c("percent.mito"))
```

'''.format(mito_gene=mito_gene, label_underline=label_underline, 
        label_stub=label_stub, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, matrix_file=matrix_file, pheno_file=pheno_file, 
        log2_transform=log2_transform, run_recluster=run_recluster, 
        run_original=run_original, descrip=descrip, cluster_type=cluster_type, 
        maximum_gene=maximum_gene, max_mito_percent=max_mito_percent)

    print '''


识别在样品之间变化差异显著的基因用于样品聚类分析 (Figure \@ref(fig:{label_stub}-mean-var-genes-seurat-plot))。首先计算每个基因的表达均值和离散度（标准差）；然后把基因按照表达值分为20个区间，计算每个区间内方差的Z-score，最小化低表达基因的检测波动引入的基因表达值的变化。

```{{r}}
#seurat_ct_{label_underline}_scale <- seurat_ct_{label_underline}@scale.data
#seurat_ct_{label_underline}_scale <- as(seurat_ct_{label_underline}_scale, matrix)
#write.table(seurat_ct_{label_underline}_scale, file=paste0(output_dir,prefix,".{label_underline}.seurat_{cluster_type}_scale_expr.xls"), sep="\\t", row.names=T, col.names=T, quote=F)
```

```{{r}}
y.cutoff_{label_underline} <- 0.5
x.low.cutoff_{label_underline} <- 0.6
x.high.cutoff_{label_underline} <- 3
```

```{{r {label_stub}-mean-var-genes-seurat-compute, include=FALSE, fig.cap="Mean-variability plot for all genes. Highly variable genes were labelled with their names."}}
seurat_ct_{label_underline} <- Seurat::MeanVarPlot(seurat_ct_{label_underline},y.cutoff = y.cutoff_{label_underline}, x.low.cutoff = x.low.cutoff_{label_underline}, fxn.x = Seurat::expMean,fxn.y = Seurat::logVarDivMean, x.high.cutoff = x.high.cutoff_{label_underline}, do.spike = ercc)
#seurat_var_genes <- seurat_ct_{label_underline}@var.genes
```

(ref:{label_stub}-mean-var-genes-seurat-plot) Mean-variability plot for all genes. Highly variable genes were labelled with their names. The X-axis represents the mean expression level, and for Y-axis represents the log(Variance/mean). Identifies genes that are outliers on a 'mean variability plot'. First, uses a function to calculate average expression (fxn.x) and dispersion (fxn.y) for each gene. Next, divides genes into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin. The purpose of this is to identify variable genes while controlling for the strong relationship between variability and average expression. 

```{{r {label_stub}-mean-var-genes-seurat-plot, fig.cap="(ref:{label_stub}-mean-var-genes-seurat-plot)"}}
seurat_ct_{label_underline} <- Seurat::MeanVarPlot(seurat_ct_{label_underline},y.cutoff = y.cutoff_{label_underline}, x.low.cutoff = x.low.cutoff_{label_underline}, fxn.x = Seurat::expMean,fxn.y = Seurat::logVarDivMean, x.high.cutoff = x.high.cutoff_{label_underline}, do.spike = ercc, do.recalc = FALSE)
seurat_var_genes_{label_underline} <- seurat_ct_{label_underline}@var.genes
```

共鉴定出`r I(length(seurat_var_genes_{label_underline}))`个样品间差异大的基因，用于区分细胞的基因表达异质性，选取其中四个绘制如 Figure \@ref(fig:{label_stub}-mean-var-genes-seurat-4)。该步筛选表达差异大的基因时未考虑样品的分组信息，部分差异大的基因可能是“组内”差异大，也有可能组间差异大，对于样品的再分类会很有帮助。

```{{r {label_stub}-mean-var-genes-seurat-4, fig.height=16, fig.cap="Violin-plot for 4 highly variable genes. Y-axis represents log-transformed expression value.", eval=T}}
#print(paste0("There are ", length(seurat_var_genes_{label_underline}), "highly variable genes identified."))
seurat_var_genes_{label_underline}_plot <- seurat_var_genes_{label_underline}[1:4]
seurat_var_genes_{label_underline}_plot <- seurat_var_genes_{label_underline}_plot[!is.na(seurat_var_genes_{label_underline}_plot)]
if(length(seurat_var_genes_{label_underline}_plot) > 1){{
    Seurat::VlnPlot(seurat_ct_{label_underline}, seurat_var_genes_{label_underline}_plot, size.x.use = 5, size.y.use=9, size.title.use = 12, size.use = 0.8, nCol=1)
}}
```

利用PCA方法根据高表达变化基因对样品进行聚类分组，结果如Figure \@ref(fig:{label_stub}-pca-plot-based-on-seurat-var-genes)所示。

```{{r {label_stub}-pca-plot-based-on-seurat-var-genes, fig.height=7, fig.cap="PCA plot of sample correlation based on the most variable genes. Each point represents one sample. Different point colors represent different sample groups. Different point shape represents different sequencing batches."}}
# Default using object@var.genes for PCA input
seurat_ct_{label_underline} <- Seurat::PCA(seurat_ct_{label_underline}, do.print=F)
pca_dim_first_{label_underline} <- dim(seurat_ct_{label_underline}@pca.x)[2]
Seurat::PCAPlot(seurat_ct_{label_underline},1,2,pt.size=2, pt.shape=pheno_data_col[2], do.return=TRUE) + scale_shape_manual(values=shapes)
```

Figure \@ref(fig:{label_stub}-seurat-vizpca) 展示了与主成分1 (PC1)和主成分2 (PC2) 最相关的前40个基因。横轴代表基因在对应主成分的rank得分，纵轴为每个基因。在PC1中得分为正的基因，其表达值在Figure \@ref(fig:{label_stub}-pca-plot-based-on-seurat-var-genes) 中右侧的样品中表达高；在PC1中得分为负的基因，其表达值在Figure \@ref(fig:{label_stub}-pca-plot-based-on-seurat-var-genes) 中左侧的样品中表达高；在PC2中得分为正的基因，其表达值在Figure \@ref(fig:{label_stub}-pca-plot-based-on-seurat-var-genes) 中上部分的样品中表达高；在PC2中得分为负的基因，其表达值在Figure \@ref(fig:{label_stub}-pca-plot-based-on-seurat-var-genes) 中下部分的样品中表达高；


```{{r {label_stub}-seurat-vizpca, fig.height=7, fig.width=7, fig.cap="Visualize top 40 genes associated with principal components. These genes are thought as classification markers for samples. Genes positively contributed to principal component 1 (PC1) highly expressed in samples at the right part of PCA plot. Genes positively controbuted to principal component2 (PC2) highly expressed in samples at the top part of PCA plot."}}
Seurat::VizPCA(seurat_ct_{label_underline},1:2, num.genes=40, font.size=0.7)
```

Figure \@ref(fig:{label_stub}-heatmap-pc1-seurat) 展示了6个主成分的关键基因在样品中的分布。每一列代表一个样品，横轴代表一个基因。样品和基因都根据其主成分得分排序，不同的子图的样品顺序不一定一致。颜色从红色-到黑色-到黄色表示基因的相对表达量由低到高。


```{{r {label_stub}-heatmap-pc1-seurat, fig.height=7, fig.cap="Heatmap showing sources of heterogeneity for PC1-6."}}
Seurat::PCHeatmap(seurat_ct_{label_underline},pc.use=1:6, do.balanced=T, label.columns = FALSE, remove.key = FALSE)
```

```{{r {label_stub}-Determine-statistically-significant-principal-components, eval=F}}
# Do 200 random samplings to find significant genes, each time randomly permute 1% of genes
# This returns a 'p-value' for each gene in each PC, based on how likely the gene/PC score woud have been observed by chance
# Note that in this case we get the same result with 200 or 1000 samplings, so we do 200 here for expediency
pca_dim_first_{label_underline} <- dim(seurat_ct_{label_underline}@pca.x)[2]
if (pca_dim_first_{label_underline} > 20) {{
    pca_dim_first_{label_underline} <- 20
}}
seurat_ct_{label_underline}=Seurat::JackStraw(seurat_ct_{label_underline},num.replicate = 200,do.print = FALSE)
```

Figure \@ref(fig:{label_stub}-Determine-statistically-significant-principal-components-plot) 展示了主成分贡献的显著性。贡献显著的主成分中更多的基因具有较低的p-value (即处于折线下方)。

```{{r {label_stub}-Determine-statistically-significant-principal-components-plot, fig.cap="Statistically significant principle components.", eval=F }}
# The jackStraw plot compares the distribution of P-values for each PC with a uniform distribution (dashed line)
# 'Significant' PCs will have a strong enrichment of genes with low p-values (solid curve above dashed line)
# In this case PC1-9 are strongly significant
Seurat::JackStrawPlot(seurat_ct_{label_underline},PCs = 1:pca_dim_first_{label_underline})
```

```{{r}}
project_pca_{label_underline} <- FALSE
```

```{{r {label_stub}-grow-gene-list-pca-seurat, fig.height=7, eval=project_pca_{label_underline}, include=FALSE}}
# Previous analysis was performed on < 400 variable genes. To identify a larger gene set that may drive
# biological differences, but did not pass our mean/variability thresholds, we first calculate PCA scores
# for all genes (PCA projection)
seurat_ct_{label_underline} <- Seurat::ProjectPCA(seurat_ct_{label_underline}, do.print=F, do.center=T)
# For test
# Seurat::PrintPCA(seurat_ct_{label_underline},1)
# Seurat::PCHeatmap(seurat_ct_{label_underline},pc.use=1, use.full=TRUE, do.balanced=T, remove.key = T, num.genes=40)
```

```{{r, eval=project_pca_{label_underline} }}
pca_dim_first_{label_underline} <- dim(seurat_ct_{label_underline}@pca.x)[2]
if (pca_dim_first_{label_underline} > 20) {{
    pca_dim_first_{label_underline} <- 20
}}
```

```{{r {label_stub}-old-seirat-sig-genes, eval=project_pca_{label_underline} }}
# Choose significant genes for PC1-9, allow each PC to contribute a max of 200 genes (to avoid one PC swamping the analysis)
seurat_ct_{label_underline}_sig_genes = Seurat::PCASigGenes(seurat_ct_{label_underline}, pcs.use=1:pca_dim_first_{label_underline}, pval.cut=0.01, max.per.pc = 200)
length(seurat_ct_{label_underline}_sig_genes)
```

```{{r, eval=project_pca_{label_underline} }}
# Now redo the PCA analysis with the new gene list
seurat_ct_{label_underline} <- Seurat::PCA(seurat_ct_{label_underline}, pc.genes=seurat_ct_{label_underline}_sig_genes, do.print=FALSE)

pca_dim_first_{label_underline} <- dim(seurat_ct_{label_underline}@pca.x)[2]
if (pca_dim_first_{label_underline} > 20) {{
    pca_dim_first_{label_underline} <- 20
}}
# Redo random sampling, PCs 1-11 are significant (additional PCs are borderline, but signal disappears with 1k replicates, though we do 200 here for expediency)
seurat_ct_{label_underline} =Seurat::JackStraw(seurat_ct_{label_underline},num.replicate = 200,do.print = FALSE)
```

```{{r {label_stub}-jackstraw-new-significant, fig.cap="Statistically significant principle components.", eval=project_pca_{label_underline} }}
Seurat::JackStrawPlot(seurat_ct_{label_underline}, PCs = 1:pca_dim_first_{label_underline})
```

```{{asis, echo=FALSE}}
PC selection – identifying the true dimensionality of a dataset – is an important step for Seurat, but can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-10 as a cutoff. We followed the jackStraw here, admittedly buoyed by seeing the PCHeatmap returning interpretable signals (including canonical dendritic cell markers) throughout these PCs. Though the results are only subtly affected by small shifts in this cutoff (you can test below), we strongly suggest always explore the PCs they choose to include downstream.
```



Figure \@ref(fig:{label_stub}-stanrarddeviation-pc) 展示了不同主成分对样品差异的贡献，便于筛选主成分用于后续分析。横轴代表识别出的主成分，1、2、3...。纵轴代表主成分反应的样品差异的度量，数值越大，对应的主成分反映的样品差异越多。

```{{r {label_stub}-stanrarddeviation-pc, fig.cap="Variations contributed by each principal component. "}}
Seurat::PCElbowPlot(seurat_ct_{label_underline})
```

```{{asis, echo={run_original} }}
Figure \@ref(fig:{label_stub}-seurat-tsne-plot-original) 展示了根据样品表达图谱进行的样品分类。
```


```{{r {label_stub}-seurat-tsne-plot-original, fig.height=16, fig.width=8, fig.cap="tSNE map of original data. Multiple iterations (1000) were performed to make sure the stability of clustering.", eval={run_original} }}
pca_dim_{label_underline} <- dim(seurat_ct_{label_underline}@pca.x)[2]
if (pca_dim_{label_underline} > 20) {{
    pca_dim_{label_underline} <- 20
}}
#pca_dim_{label_underline} <- 20
seurat_ct_{label_underline} <- Seurat::RunTSNE(seurat_ct_{label_underline}, dims.use=1:pca_dim_{label_underline}, k.seed=11521, do.fast=T, perplexity = round(seurat_samplenumber_{label_underline}/5))
Seurat::TSNEPlot(seurat_ct_{label_underline}, pt.size=1, do.return=T, do.label=T, no.legend=T) +
    labs(title=paste0("perplexity = ",round(seurat_samplenumber_{label_underline}/5))) +
    theme(legend.text=element_text(size=5))
```            


```{{asis, echo={run_recluster} }}
Figure \@ref(fig:{label_stub}-seurat-findclusters) 展示了根据样品表达图谱进行的不同精细程度的分类。resolution越大，分出的类越多。可以直接选择大的resolution分成小类，也可以先分出大类，然后再细分亚类。
```

```{{r {label_stub}-seurat-findclusters, fig.height=13, fig.cap="Test different parameters for density clustering to make sure get the best re-classification.", eval={run_recluster} }}
findc_resolution_{label_underline} <- data.frame(res=c(0.4,0.7,1.0))
pca_dim_{label_underline} <- dim(seurat_ct_{label_underline}@pca.x)[2]

find_cluster_tsne_list_{label_underline} <- by(data=findc_resolution_{label_underline}, INDICES = findc_resolution_{label_underline}$res, FUN=function(x){{
  res <- x$res
  seurat_ct_{label_underline} <- Seurat::FindClusters(seurat_ct_{label_underline}, pc.use=1:pca_dim_{label_underline}, resolution=res, print.output=0, random.seed=res*100, save.SNN=T, temp.file.location="./")
  ident_{label_underline}_cluster <- length(unique(Seurat::FetchData(seurat_ct_{label_underline}, c("ident"))$ident))
  if (ident_{label_underline}_cluster<4) {{
    ident_{label_underline}_cluster = 4
  }}
  perplexity = round(seurat_samplenumber_{label_underline}/ident_{label_underline}_cluster)
  if (perplexity>30){{
    perplexity = 30
  }}
  seurat_ct_{label_underline} <- Seurat::RunTSNE(seurat_ct_{label_underline}, dims.use=1:pca_dim_{label_underline}, k.seed=res*100, do.fast=T, perplexity = perplexity)
  Seurat::TSNEPlot(seurat_ct_{label_underline}, pt.size=1, do.return=T, do.label=T, no.legend=T) +
    labs(title=paste0("perplexity = ", perplexity, "; resolution = ", res)) +
    theme(legend.text=element_text(size=5))
}})


do.call(grid.arrange, c(find_cluster_tsne_list_{label_underline}, ncol=2))
```

```{{r}}
ident_{label_underline} <- unique(Seurat::FetchData(seurat_ct_{label_underline}, c("ident"))$ident)
```

```{{r {label_stub}-seurat-findclusters-final, fig.cap="Sample re-classification.", eval={run_recluster} }}
res_{label_underline}=1
seed_{label_underline}=res_{label_underline} * 100
seurat_ct_{label_underline} <- Seurat::FindClusters(seurat_ct_{label_underline}, pc.use=1:pca_dim_{label_underline}, resolution=res_{label_underline}, print.output=0, random.seed=seed_{label_underline}, save.SNN=T, temp.file.location="./")
ident_{label_underline}_cluster <- length(unique(Seurat::FetchData(seurat_ct_{label_underline}, c("ident"))$ident))
if (ident_{label_underline}_cluster<4) {{
  ident_{label_underline}_cluster = 4
}}
perplexity = round(seurat_samplenumber_{label_underline}/ident_{label_underline}_cluster)
if (perplexity>30){{
  perplexity = 30
}}
seurat_ct_{label_underline} <- Seurat::RunTSNE(seurat_ct_{label_underline}, dims.use=1:pca_dim_{label_underline}, k.seed=seed_{label_underline}, do.fast=T, perplexity = perplexity)
Seurat::TSNEPlot(seurat_ct_{label_underline}, pt.size=1, do.return=T, do.label=T, no.legend=TRUE) +
    labs(title=paste0("perplexity = ",perplexity, "; resolution = ", res_{label_underline})) +
    theme(legend.text=element_text(size=5))
seurat_identity_{label_underline} = Seurat::FetchData(seurat_ct_{label_underline}, c("ident", "orig.ident"))
ident_{label_underline} <- unique(seurat_identity_{label_underline}$ident)
write.table(seurat_identity_{label_underline}, file=paste0(output_dir,prefix,".{label_underline}.seurat_{cluster_type}_phenoData.xls"), sep="\\t", row.names=T, col.names=T, quote=F)
```


选取**resolution=`r res_{label_underline}` (random seed=`seed_{label_underline}`)**作为分类标准，进行下游探索 (Figure \@ref(fig:{label_stub}-seurat-findclusters-final))。分类结果以及新分类与原分类的对应表格可[点击下载](`r  I(paste0(link_dir, prefix,".{label_underline}.seurat_{cluster_type}_phenoData.xls"))`)。


```{{r {label_stub}-save-seurat-{cluster_type}, eval={run_recluster} }}
save(seurat_ct_{label_underline},file=paste0(prefix, ".{label_underline}.seurat_{cluster_type}.Robj"))
```

```{{r, eval={run_recluster} }}
# output ident and tsne
seurat_ct_{label_underline}_tsne_identity_pca <- cbind(ident=seurat_ct_{label_underline}@ident, seurat_ct_{label_underline}@tsne.rot, seurat_ct_{label_underline}@pca.rot[,1:4])
write.table(seurat_ct_{label_underline}_tsne_identity_pca, file=paste0(output_dir,prefix,".{label_underline}.seurat_{cluster_type}_coord.xls"), sep="\\t", row.names=T, col.names=T, quote=F)
```

```{{r, eval={run_recluster} }}
BuildClusterTree(seurat_ct_{label_underline}, pcs.use=1:pca_dim_{label_underline})
```

```{{r}}
findMarker_{label_underline} <- FALSE
```

```{{r {label_stub}-find-all-markers, eval=findMarker_{label_underline} }}
# thresh.use: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. 
# min.pct: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. 
all_seurat_markers_{label_underline} = Seurat::FindAllMarkers(seurat_ct_{label_underline}, thresh.use=0.25, min.pct=0.2, test.use="roc", print.bar = debug, only.pos=TRUE)
```

```{{asis, echo=findMarker_{label_underline} }}
样品标记基因列表 (Table \@ref(tab:{label_stub}-list-marker-genes))。

Table: (\#tab:{label_stub}-marker-gene-label-explanation) Explanation for each column of marker files.
```

```{{r {label_stub}-marker-gene-label-explanation, results="asis", eval=findMarker_{label_underline} }}
marker_gene_label_explanation <- "Name;Explanation
myAUC;Area under ROC (Receiver Operating Characteristic) curve, the larger the better
avg_diff;X-fold difference (log-scale)
power;Classification power for each gene, ranging from 0 (random) - 1 (perfect). Though not a statistical test, it is often very useful for finding clean markers.
pct.1;Percentage of samples have this gene detected in target sample class
pct.2;Percentage of samples have this gene detected in other sample classes
cluster;Target sample class
gene;Gene name
"
marker_gene_label_explanation_data <- read.table(text=marker_gene_label_explanation, sep=";",header=T,quote="")
knitr::kable(marker_gene_label_explanation_data, booktabs=T, format="markdown")
```

```{{asis, echo=findMarker_{label_underline} }}
**样品标记基因示例**
```

```{{r {label_stub}-sample-marker-gene-example-top5, eval=findMarker_{label_underline} }}
sample_marker_gene_exp_top5_{label_underline} <- all_seurat_markers_{label_underline} %>% group_by(cluster) %>% top_n(5, avg_diff)
knitr::kable(sample_marker_gene_exp_top5_{label_underline}, booktabs=T, caption="List of top 5 markers for each cluster.")
```

```{{asis, echo=findMarker_{label_underline} }}
Table: (\#tab:{label_stub}-list-marker-genes) List of markers for each sample class. Currently only positive genes are reported to shrink computation time.
```

```{{r {label_stub}-list-marker-genes, results="asis", eval=findMarker_{label_underline} }}
ident2_{label_underline} <- ident_{label_underline}[order(ident_{label_underline})]

listMarkers_{label_underline} <- function(cluster){{
  c(paste0("[Cluster ", cluster," UP marker](",paste0(link_dir, prefix,".{label_underline}.seurat.",cluster,'.up_marker.xls'),")"),
  paste0("[Cluster ", cluster," DW marker](",paste0(link_dir, prefix,".{label_underline}.seurat.",cluster,'.dw_marker.xls'),")"))
}}
marker_list_tmp_{label_underline} <- as.data.frame(t(sapply(ident2_{label_underline}, listMarkers_{label_underline})))
colnames(marker_list_tmp_{label_underline}) <- c("UP","DW (ignore)")
#rownames(marker_list_tmp_{label_underline}) <- ident2_{label_underline}
knitr::kable(marker_list_tmp_{label_underline}, booktabs=T, format="markdown", caption="List of markers for each sample class.")
```


```{{r {label_stub}-output-all-markers, eval=findMarker_{label_underline} }}
outputMarkers_{label_underline} <- function(cluster){{
  up_marker <- all_seurat_markers_{label_underline}[all_seurat_markers_{label_underline}$cluster==cluster & all_seurat_markers_{label_underline}$avg_diff>0,]
  write.table(up_marker, file=paste0(output_dir,prefix,".{label_underline}.seurat.",cluster,'.up_marker.xls'), sep="\\t", row.names=F, col.names=T, quote=F)
  dw_marker <- all_seurat_markers_{label_underline}[all_seurat_markers_{label_underline}$cluster==cluster & all_seurat_markers_{label_underline}$avg_diff<0,]
  write.table(dw_marker, file=paste0(output_dir,prefix,".{label_underline}.seurat.",cluster,'.dw_marker.xls'), sep="\\t", row.names=F, col.names=T, quote=F)
  c(up_marker$gene[1], up_marker$gene[2])
}}
test_marker_for_plot_{label_underline} <- sapply(ident_{label_underline}, outputMarkers_{label_underline})
```

```{{asis, echo=findMarker_{label_underline} }}
Figure \@ref(fig:{label_stub}-plot-representative-markers) 展示了代表性样品Marker基因在对应样品的表达图谱。
```

```{{r {label_stub}-plot-representative-markers, fig.height=18, fig.cap="Violin-plot for representative markers. Y-axis represents log-transformed expression value.", eval=findMarker_{label_underline} }}
test_marker_for_plot_{label_underline} <- test_marker_for_plot_{label_underline}[!is.na(test_marker_for_plot_{label_underline})]
Seurat::VlnPlot(seurat_ct_{label_underline}, test_marker_for_plot_{label_underline}, size.x.use = 5, size.y.use=6, size.title.use = 10, size.use = 0.8, nCol=2)
```

```{{asis, echo=findMarker_{label_underline} }}
Figure \@ref(fig:{label_stub}-seurat-featureplot) 展示数个样品标记基因在样品间的表达情况（灰色到蓝色表达由低到高）。
```

```{{r {label_stub}-seurat-featureplot, fig.height=7, fig.cap="Map of marker genes to related sample classifications. At most two marker genes were selected for each sample.  Different colors represent the expression (low:grey;high:blue) of the gene in the context of all the cells and help validate the specificity of the marker.", eval=findMarker_{label_underline} }}
Seurat::FeaturePlot(seurat_ct_{label_underline}, test_marker_for_plot_{label_underline},cols.use = c("grey","dark blue"))
```

```{{asis, echo=findMarker_{label_underline} }}
Figure \@ref(fig:{label_stub}-seurat-heatmap-20-markers) 展示了每个cluster最具代表性的10个基因的相对表达图谱 （颜色从红色-黑色-黄色代表相对表达值由低到高）。
```

```{{r {label_stub}-seurat-heatmap-20-markers, fig.height=14, fig.cap="Heatmap showing expression patterns of top 20 marker genes for each group. Colors from read to black to yellow showing Z-scaled expression value from low to high.", eval=findMarker_{label_underline} }}
all_seurat_markers_{label_underline} %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_diff) -> top10_{label_underline}
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
Seurat::DoHeatmap(seurat_ct_{label_underline}, genes.use = top10_{label_underline}$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)
```

```{{r}}
system(paste0("mkdir -p ", link_dir))
system(paste0("/bin/cp -up ", output_dir,"* ", link_dir))
```

'''.format(label=label, label_underline=label_underline, label_stub=label_stub, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, matrix_file=matrix_file, pheno_file=pheno_file, 
        log2_transform=log2_transform, run_recluster=run_recluster, 
        run_original=run_original, descrip=descrip, cluster_type=cluster_type)
    #-----------end close fh-----------

#```{{r {label_stub}-seurat-tsne-plot-original-test, fig.height=16, fig.width=8, fig.cap="tSNE map of original data. Multiple iterations (1000) were performed to make sure the stability of clustering.", eval=FALSE }}
#pca_dim_{label_underline} <- dim(seurat_ct_{label_underline}@pca.x)[2]
#plot_tsne_list_{label_underline} <- by(data=cluster_perplexity_{label_underline}, INDICES = cluster_perplexity_{label_underline}$a, FUN= function(x) {{
#  x <- x$a
#  seurat_ct_{label_underline} <- Seurat::RunTSNE(seurat_ct_{label_underline}, dims.use=1:pca_dim_{label_underline}, k.seed=11521, max_iter=1000, perplexity = x, do.fast=T)
#  Seurat::TSNEPlot(seurat_ct_{label_underline}, pt.size=1, do.return=T, no.legend=TRUE, do.label=T) +
#    labs(title=paste0("perplexity = ",x)) + theme(legend.text=element_text(size=5))
#}})
#
#do.call(grid.arrange, c(plot_tsne_list_{label_underline}, ncol=2))
#```

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


