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
        help="Filtered RDS file. Outlier samples and genes removed, only endog_genes kept. <exprs> slot containing final expression data.")
    parser.add_option("-l", "--label", dest="label",
        help="several concise words to represent the content of this analysis. Only alphabets and blanks are allowed. Must be unique if this script will be run multiple times. String like <scran norm original cluster>, <cpm norm reclassify>, <scran norm batch remove reclassify>,  <cpm norm sc3 cluster>")
    parser.add_option("-m", "--matrix-file", dest="matrix_file",
        help="A matrix file containing **normalized** **log2 transformed** gene expression according to format specified above. One of <-r> and <-m> is required. If both supplied, <-m> would have higher priority.")
    parser.add_option("-p", "--pheno-file", dest="pheno_file",
        help="(lowercase p) A pheno file according to format specified above. OPTIONAL.")
    parser.add_option("-d", "--description", dest="descrip",
        help="A paragraph of markdown formated words to describe things one wants to do.")
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
    assert options.label != None, "A label needed for -l"
    assert options.descrip != None, "A label needed for -d"
    assert (options.rds_file != None) or (options.matrix_file != None), "A rds or matrix file needed for -r or -m"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    label = options.label
    label_underline = label.replace(' ', '_')
    label_stub = label.replace(' ', '-') 
    output_directory = options.output_directory
    link_directory   = options.link_directory
    prefix   = options.prefix
    matrix_file = options.matrix_file
    pheno_file  = options.pheno_file
    rds_file = options.rds_file
    descrip  = options.descrip
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
    print '''
```{{r, include=F}}

if(! exists("output_dir")){{
    output_dir <- "{output_directory}"
}}

if(! exists("link_dir")){{
    link_dir <- "{link_directory}"
}}

if(! exists("prefix")){{
    prefix <- "{prefix}"
}}

library(scater, quietly = TRUE)

```

# {descrip} ({label})


{descrip}
'''.format(label=label, label_underline=label_underline, label_stub=label_stub, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, descrip=descrip)

    if matrix_file:
        print '''
```{{r initiate-sc3-{label_stub}-matrix, cache.extra=mtime("{matrix_file}")}}
matrix_data <- read.table("{matrix_file}", sep="\\t", row.names=1, header=T, quote="", check.names=F)
```

'''.format(label_stub=label_stub, matrix_file=matrix_file, pheno_file=pheno_file)
        if pheno_file:
            print '''
```{{r initiate-sc3-{label_stub}-matrix-pheno, cache.extra=mtime("{pheno_file}")}}
pheno_data <- read.table("{pheno_file}", sep="\\t", header=T, row.names=1, quote="", check.names=F)
pheno_data <- pheno_data[match(colnames(matrix_data), rownames(pheno_data)), ]
pheno_data_col = colnames(pheno_data)
pheno_data_col_len <- length(pheno_data_col)
if(pheno_data_col_len<2) {{
  pheno_data_col <- rep(pheno_data_col, 2)
}}
  
if((pheno_data_col[1]!="conditions") | (pheno_data_col[2]!="batch")) {{
  stop("***Please check if phenoData has conditions and batch at the second and third column.***")
}}

pheno_dataframe <- new("AnnotatedDataFrame",  pheno_data)

sceset_data_{label_underline}_qc_sc3 <- scater::newSCESet(
  exprsData = matrix_data, 
  phenoData = pheno_dataframe
  )
sceset_data_{label_underline}_qc_sc3@logged = TRUE
```

'''.format(label_stub=label_stub, matrix_file=matrix_file, pheno_file=pheno_file, 
        label_underline=label_underline)
        else:
            print '''
```{{r initiate-sc3-{label_stub}-matrix-2 }}
sceset_data_{label_underline}_qc_sc3 <- scater::newSCESet(
  exprsData = matrix_data
  )
sceset_data_{label_underline}_qc_sc3@logged = TRUE
```
'''.format(label_stub=label_stub, label_underline=label_underline)

    else:
        print '''

```{{r initiate-sc3-{label_stub}-data, cache.extra=mtime("{rds_file}")}}
sceset_data_{label_underline}_qc_sc3 <- readRDS("{rds_file}")
featureNames(sceset_data_{label_underline}_qc_sc3) <- featureData(sceset_data_{label_underline}_qc_sc3)$Associated_Gene_Name
```
'''.format(label_stub=label_stub, label_underline=label_underline, rds_file=rds_file)

    print '''
```{{r initiate-sc3-{label_stub}-attr}}
a_{label_underline} = dim(sceset_data_{label_underline}_qc_sc3)
geneNumber_{label_underline} = as.vector(a_{label_underline}[1]) 
sampleNumber_{label_underline} = as.vector(a_{label_underline}[2])


## Uncomment when necessary
# level <- length(unique(unlist(pData(sceset_data_{label_underline}_qc_sc3)[pheno_data_col[2]])))
# shapes = (1:level)%%30  # maximum allow 30 types of symbols
# add following line to ggplot2 plot
#+ scale_shape_manual(values=shapes)
```

```{{r sc3-{label_stub}, include=FALSE }}
sceset_data_{label_underline}_qc_sc3 <- SC3::sc3_prepare(sceset_data_{label_underline}_qc_sc3, exprs_values="exprs", ks=3, gene_filter = F, n_cores = 30, svm_max=3000, rand_seed=11521)
sceset_data_{label_underline}_qc_sc3 <- SC3::sc3_estimate_k(sceset_data_{label_underline}_qc_sc3)

sampleClusterNumber_{label_underline} <- sceset_data_{label_underline}_qc_sc3@sc3$k_estimation
sceset_data_{label_underline}_qc_sc3 <- SC3::sc3(sceset_data_{label_underline}_qc_sc3, exprs_values="exprs", ks=sampleClusterNumber_{label_underline}, gene_filter = F, biology=T, n_cores = 30, k_estimator=F, svm_max=3000, rand_seed=11521)
sceset_data_{label_underline}_qc_sc3_result_suf = ".sc3_{label_underline}_cluster_marker_gene_total_results.xls"
SC3::sc3_export_results_xls(sceset_data_{label_underline}_qc_sc3,
      filename=paste0(output_dir,prefix, sceset_data_{label_underline}_qc_sc3_result_suf))
```

```{{r eval=debug}}
# Interactively select the best cluster number
#SC3::sc3_interactive(sceset_data_{label_underline}_qc_sc3)
```

使用软件SC3进行聚类分析，数据评估计算的最优分类数为`r I(sampleClusterNumber_{label_underline})`，聚类结果如 Figure \@ref(fig:sc3-{label_stub}-pca) and \@ref(fig:sc3-{label_stub}-tsne)。总计算结果、差异基因列表、Marker基因列表存储于 [`r I(paste0(prefix,sceset_data_{label_underline}_qc_sc3_result_suf))`](`r I(paste0(link_dir, prefix, sceset_data_{label_underline}_qc_sc3_result_suf))`).


```{{r sc3-{label_stub}-plot}}
sceset_data_{label_underline}_qc_sc3_cluster_colname <- paste0("sc3_",sampleClusterNumber_{label_underline},"_clusters")
pData(sceset_data_{label_underline}_qc_sc3)[sceset_data_{label_underline}_qc_sc3_cluster_colname] <- lapply(pData(sceset_data_{label_underline}_qc_sc3)[sceset_data_{label_underline}_qc_sc3_cluster_colname], as.factor)
```

```{{r sc3-{label_stub}-txt-data}}
sceset_data_{label_underline}_qc_sc3_pData <- pData(sceset_data_{label_underline}_qc_sc3)
sc3_{label_underline}_cluster = sceset_data_{label_underline}_qc_sc3_pData[,grep(paste0("sc3_", sampleClusterNumber_{label_underline}), colnames(sceset_data_{label_underline}_qc_sc3_pData))]
write.table(sc3_{label_underline}_cluster, file=paste0(output_dir,prefix,".sc3_{label_underline}_cluster_results.xls"), sep="\\t", row.names = T, col.names=T, quote=F)

sceset_data_{label_underline}_qc_sc3_fData <- fData(sceset_data_{label_underline}_qc_sc3)
sceset_data_{label_underline}_qc_sc3_marker = sceset_data_{label_underline}_qc_sc3_fData[,grep(paste0("sc3_",sampleClusterNumber_{label_underline}), colnames(sceset_data_{label_underline}_qc_sc3_fData))]
write.table(sceset_data_{label_underline}_qc_sc3_marker, file=paste0(output_dir,prefix,".sc3_{label_underline}_marker_gene_results.xls"), sep="\\t", row.names = T, col.names=T, quote=F)

sceset_data_{label_underline}_qc_sc3_expr <- exprs(sceset_data_{label_underline}_qc_sc3)
sceset_data_{label_underline}_qc_sc3_expr_colname <- colnames(sceset_data_{label_underline}_qc_sc3_expr)
sceset_data_{label_underline}_qc_sc3_expr_cluster <- sc3_{label_underline}_cluster[,1]
colnames(sceset_data_{label_underline}_qc_sc3_expr) <- paste(sceset_data_{label_underline}_qc_sc3_expr_colname,sceset_data_{label_underline}_qc_sc3_expr_cluster,sep="___")
write.table(sceset_data_{label_underline}_qc_sc3_expr, file=paste0(output_dir,prefix,".sc3_{label_underline}_reclustered_expr_mat.xls"), sep="\\t", row.names=T, col.names=T, quote=F)

sceset_data_{label_underline}_qc_sc3_pheno <- pData(sceset_data_{label_underline}_qc_sc3)

sceset_data_{label_underline}_qc_sc3_pheno <- cbind(sceset_data_{label_underline}_qc_sc3_pheno[,1:pheno_data_col_len], sc3_cluster=sceset_data_{label_underline}_qc_sc3_pheno[, grep(sceset_data_{label_underline}_qc_sc3_cluster_colname, colnames(sceset_data_{label_underline}_qc_sc3_pheno))])

sceset_data_{label_underline}_qc_sc3_pheno$id <- rownames(sceset_data_{label_underline}_qc_sc3_pheno)
rownames(sceset_data_{label_underline}_qc_sc3_pheno) <- colnames(sceset_data_{label_underline}_qc_sc3_expr)
write.table(sceset_data_{label_underline}_qc_sc3_pheno, file=paste0(output_dir,prefix,".sc3_{label_underline}_reclustered_pheno.xls"),sep="\\t",row.names=T,col.names=T,quote=F)



sceset_data_{label_underline}_qc_sc3_outputMarkers <- function(cluster) {{
  marker <- sceset_data_{label_underline}_qc_sc3_marker[ sceset_data_{label_underline}_qc_sc3_marker[paste0("sc3_",sampleClusterNumber_{label_underline},"_markers_clusts")]==cluster & 
                                  sceset_data_{label_underline}_qc_sc3_marker[paste0("sc3_",sampleClusterNumber_{label_underline},"_markers_padj")]<=0.01 &
                                  sceset_data_{label_underline}_qc_sc3_marker[paste0("sc3_",sampleClusterNumber_{label_underline},"_markers_auroc")]<=0.85,]
  write.table(marker, file=paste0(output_dir,prefix,".sc3_{label_underline}_cluster_", cluster, "_markers.xls"), sep="\\t", row.names=T, col.names=T, quote=F)
}}

sceset_data_{label_underline}_qc_sc3_outputMarkers_tmp <- sapply(1:sampleClusterNumber_{label_underline}, sceset_data_{label_underline}_qc_sc3_outputMarkers)
```

Table: (\#tab:sc3-{label_stub}-marker-gene) List of markers for each sample class (padj<=0.01 && auroc<=0.85).

```{{r sc3-{label_stub}-marker-gene}}
sceset_data_{label_underline}_qc_sc3_listMarkers <- function(cluster) {{
  c(paste0("[Cluster ", cluster," marker](",paste0(link_dir, prefix,".sc3_{label_underline}_cluster_", cluster, "_markers.xls"),")"))
}}
sceset_data_{label_underline}_qc_sc3_listMarkers_df <- as.data.frame(sapply(1:sampleClusterNumber_{label_underline}, sceset_data_{label_underline}_qc_sc3_listMarkers))
colnames(sceset_data_{label_underline}_qc_sc3_listMarkers_df) <- c("Markers")
#rownames(sceset_data_{label_underline}_qc_sc3_listMarkers_df) <- 1:sampleClusterNumber_{label_underline}
knitr::kable(sceset_data_{label_underline}_qc_sc3_listMarkers_df, booktabs=T, format="markdown", caption="List of markers for each sample class.")
```

```{{r}}
perplexity_{label_underline} <- sampleClusterNumber_{label_underline}
if (perplexity_{label_underline}<5){{
    perplexity_{label_underline} = 5
}}
perplexity_{label_underline} = round(sampleNumber_{label_underline}/perplexity_{label_underline})
```

```{{r sc3-{label_stub}-pca, fig.cap="PCA map showing sample classification with each sample colored by assigned new clusters."}}
scater::plotPCA(sceset_data_{label_underline}_qc_sc3, 
                ntop=5000, 
                colour_by=sceset_data_{label_underline}_qc_sc3_cluster_colname,
                size_by=sceset_data_{label_underline}_qc_sc3_cluster_colname,
                #size_by=paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score"),
                exprs_values="exprs") 
```

```{{r sc3-{label_stub}-pca2, fig.width=12, fig.height=12, fig.cap="PCA map showing sample classification with each sample colored by assigned new clusters."}}
scater::plotPCA(sceset_data_{label_underline}_qc_sc3, 
                ntop=5000, ncomponents=3, 
                colour_by=sceset_data_{label_underline}_qc_sc3_cluster_colname,
                size_by=sceset_data_{label_underline}_qc_sc3_cluster_colname,
                #size_by=paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score"),
                exprs_values="exprs") 
```

```{{r sc3-{label_stub}-tsne, fig.cap="tSNE map showing sample classification with each sample colored by assigned new clusters."}}
scater::plotTSNE(sceset_data_{label_underline}_qc_sc3, 
         ntop=5000, perplexity = perplexity_{label_underline},
         colour_by=sceset_data_{label_underline}_qc_sc3_cluster_colname,
         size_by=sceset_data_{label_underline}_qc_sc3_cluster_colname,
         exprs_values="exprs",
         rand_seed = 11521) + labs(title=paste0("perplexity = ",perplexity_{label_underline}))+guides(fill=guide_legend(ncol=2)) 
         #+ geom_text(aes_string(label=sceset_data_{label_underline}_qc_sc3_cluster_colname), check_overlap=T)
```

```{{r sc3-{label_stub}-tsne2, fig.cap="tSNE map showing sample classification with each sample colored by assigned new clusters."}}
scater::plotTSNE(sceset_data_{label_underline}_qc_sc3, 
         ntop=5000, perplexity = perplexity_{label_underline},
         colour_by=pheno_data_col[1],
         size_by=sceset_data_{label_underline}_qc_sc3_cluster_colname,
         exprs_values="exprs",
         rand_seed = 11521) + labs(title=paste0("perplexity = ",perplexity_{label_underline}))+guides(fill=guide_legend(ncol=2)) 
         #+ geom_text(aes_string(label=sceset_data_{label_underline}_qc_sc3_cluster_colname), check_overlap=T)
```

Figure \@ref(fig:sc3-{label_stub}-consensus) 为聚类一致性矩阵，每一行和每一列都为一个细胞，为对称矩阵。每一个方格代表对应细胞之间的相似性 （相似性的计算方式是这两个细胞在不同聚类参数下都聚为同一个类的比例）。相似性为0 （蓝色）表示这两个细胞总是被分在不同的类里面。相似性为1（红色）表示这两个细胞总是聚在一个类里面。最理想的聚类结果是所有对角框为红色，其它框为蓝色。上端的颜色条纹代表样品分组，从左到右为Cluster1, 2, 3....。

(ref:sc3-{label_stub}-consensus) The consensus matrix is a N by N matrix, where N is the number of cells in the input dataset. It represents similarity between the cells based on the averaging of clustering results from all combinations of clustering parameters. Similarity 0 (blue) means that the two cells are always assigned to different clusters. In contrast, similarity 1 (red) means that the two cells are always assigned to the same cluster. The consensus matrix is clustered by hierarchical clustering and has a diagonal-block structure. Intuitively, the perfect clustering is achieved when all diagonal blocks are completely red and all off-diagonal elements are completely blue. [PDF](`r I(paste0(link_dir, prefix, "sc3_{label_underline}_consensus.pdf"))`)

```{{r sc3-{label_stub}-consensus-plot, include=F}}
sc3_{label_underline}_consensus_pic = paste0(output_dir, prefix, "sc3_{label_underline}_consensus.")
pdf(file=paste0(sc3_{label_underline}_consensus_pic, "pdf"), width=20, height=20, pointsize=10)
SC3::sc3_plot_consensus(sceset_data_{label_underline}_qc_sc3, k=sampleClusterNumber_{label_underline}, show_pdata=c(pheno_data_col,sceset_data_{label_underline}_qc_sc3_cluster_colname, paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score")))
dev.off()
png(file=paste0(sc3_{label_underline}_consensus_pic, "png"), width=27, height=27, pointsize=9, units="cm", res=300)
SC3::sc3_plot_consensus(sceset_data_{label_underline}_qc_sc3, k=sampleClusterNumber_{label_underline}, show_pdata=c(pheno_data_col,sceset_data_{label_underline}_qc_sc3_cluster_colname, paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score")))
dev.off()
```

```{{r sc3-{label_stub}-consensus, fig.cap="(ref:sc3-{label_stub}-consensus)"}}
knitr::include_graphics(paste0(link_dir, prefix, "sc3_{label_underline}_consensus.png"))
```


Figure \@ref(fig:sc3-{label_stub}-deg) 列出了Top 50差异表达的基因。上端的颜色条纹代表样品分组，从左到右为Cluster1, 2, 3....。

(ref:sc3-{label_stub}-deg) Differential expression is calculated using the non-parametric Kruskal-Wallis test. A significant p-value indicates that gene expression in at least one cluster stochastically dominates one other cluster. SC3 provides a list of all differentially expressed genes with adjusted p-values < 0.01 and plots gene expression profiles of the 50 genes with the lowest p-values. Note that the calculation of differential expression after clustering can introduce a bias in the distribution of p-values, and thus we advise to use the p-values for ranking the genes only. [PDF](`r I(paste0(link_dir, prefix, "sc3_{label_underline}_deg.pdf"))`)

```{{r sc3-{label_stub}-deg-plot, include=F}}
pdf(file=paste0(output_dir, prefix, "sc3_{label_underline}_deg.pdf"), width=15, height=25, pointsize=10)
SC3::sc3_plot_de_genes(sceset_data_{label_underline}_qc_sc3, 
                       k=sampleClusterNumber_{label_underline}, 
                       show_pdata=c(sceset_data_{label_underline}_qc_sc3_cluster_colname,
                                    paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score")))
dev.off()
png(file=paste0(output_dir, prefix, "sc3_{label_underline}_deg.png"), width=30, height=25, pointsize=10, units="cm", res=300)
SC3::sc3_plot_de_genes(sceset_data_{label_underline}_qc_sc3, 
                       k=sampleClusterNumber_{label_underline}, 
                       show_pdata=c(sceset_data_{label_underline}_qc_sc3_cluster_colname,
                                    paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score")))
dev.off()
```

```{{r sc3-{label_stub}-deg, fig.cap="(ref:sc3-{label_stub}-deg)"}}
knitr::include_graphics(paste0(link_dir, prefix, "sc3_{label_underline}_deg.png"))
```


Figure \@ref(fig:sc3-{label_stub}-marker) 列出了每个类的top 10的标记基因。图左侧的颜色块对应于Marker基因所在的类。其顺序与图例中`legend`的标记顺序一致。

(ref:sc3-{label_stub}-marker) To find marker genes, for each gene a binary classifier is constructed based on the mean cluster expression values. The classifier prediction is then calculated using the gene expression ranks. The area under the receiver operating characteristic (ROC) curve is used to quantify the accuracy of the prediction. A p-value is assigned to each gene by using the Wilcoxon signed rank test. By default the genes with the area under the ROC curve (AUROC) > 0.85 and with the p-value < 0.01 are selected and the top 10 marker genes of each cluster are visualized in this heatmap. [PDF](`r I(paste0(link_dir, prefix, "sc3_{label_underline}_marker.pdf"))`)

```{{r sc3-{label_stub}-marker-plot, include=F}}
pdf(file=paste0(output_dir, prefix, "sc3_{label_underline}_marker.pdf"), width=17, height=25, pointsize=9)
SC3::sc3_plot_markers(sceset_data_{label_underline}_qc_sc3, 
                       k=sampleClusterNumber_{label_underline}, 
                       show_pdata=c(paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score"), pheno_data_col))
dev.off()

png(file=paste0(output_dir, prefix, "sc3_{label_underline}_marker.png"), width=35, height=50, pointsize=9, units="cm", res=300)
SC3::sc3_plot_markers(sceset_data_{label_underline}_qc_sc3, 
                       k=sampleClusterNumber_{label_underline}, 
                       show_pdata=c(paste0("sc3_",sampleClusterNumber_{label_underline},"_log2_outlier_score"), pheno_data_col))
dev.off()
```

```{{r sc3-{label_stub}-marker, fig.cap="(ref:sc3-{label_stub}-marker)"}}
knitr::include_graphics(paste0(link_dir, prefix, "sc3_{label_underline}_marker.png"))
```

```{{r}}
system(paste0("mkdir -p ", link_dir))
system(paste0("/bin/cp -up ", output_dir,"* ", link_dir))
```

'''.format(rds_file=rds_file, label=label, label_underline=label_underline, label_stub=label_stub, 
        output_directory=output_directory, link_directory=link_directory, 
        prefix=prefix, descrip=descrip)
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


