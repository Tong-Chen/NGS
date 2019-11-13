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

This script is used to do seurat analysis.

${txtbld}OPTIONS${txtrst}:
	-f	Data file ${bldred}[NECESSARY]${txtrst}
	-z	Is there a header[${bldred}Default TRUE${txtrst}]
EOF
}

file=
header='TRUE'

while getopts "hf:z:" OPTION
do
	case $OPTION in
		h)
			usage
			exit 1
			;;
		f)
			file=$OPTARG
			;;
		z)
			header=$OPTARG
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

cat <<END >

## 细胞分型和生物标志物筛选 (Seurat & normalized batch corrected data)

Seurat是纽约基因组中心的SATIJA教授实验室开发的用于单细胞转录组数据质控、分析和探索的工具包，发表于Cell、NBT等杂志。可以用来识别细胞中高度变化的基因、差异表达的基因和细胞类型标记基因。并可利用PCA和t-SNE方法鉴定细胞类型、细胞状态和细胞空间定位。

质控过滤并且标准化移除批次效应的数据导入Seurat。

```{asis,echo=FALSE}
Regress out unwanted sources of variation

Your single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. Seurat implements a basic version of this by constructing linear models to predict gene expression based on user-defined variables. Seurat stores the z-scored residuals of these models in the scale.data slot, and they are used for dimensionality reduction and clustering.

We typically regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data), the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a ‘cell-cycle’ score (as in Macosko et al) and Regress this out as well. In this simple example here for post-mitotic blood cells, we simply regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content an example.

#note that this overwrites pbmc@scale.data. Therefore, if you intend to use RegressOut, you can set do.scale=F and do.center=F in the original object to save some time.
pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
```

```{r scran-norm-batchlmR-data-log-transformed-for-seaurt}
exprs_batch_corrected_scran_colname = colnames(exprs_batch_corrected_scran)
conditions_tmp_scran_norm_batchlmR_br = pData(sceset_data_scran_qc)[[pheno_data_col[1]]]
exprs_batch_corrected_scran_colname_seurat = paste(exprs_batch_corrected_scran_colname,conditions_tmp_scran_norm_batchlmR_br,sep="@")
exprs_batch_corrected_scran_seurat <- exprs_batch_corrected_scran
colnames(exprs_batch_corrected_scran_seurat) <- exprs_batch_corrected_scran_colname_seurat
rownames(exprs_batch_corrected_scran_seurat) <- featureData(sceset_data_scran_qc)$Associated_Gene_Name
write.table(exprs_batch_corrected_scran_seurat, file=paste0(output_dir,prefix,".seurat.scran_norm_log2transform_batchlmR_batch_rmove.xls"),sep="\t",row.names=T,col.names=T,quote=F)
seurat_pheno_scran_norm_batchlmR_br = pData(sceset_data_scran_qc)[,1:pheno_data_col_len]
seurat_pheno_scran_norm_batchlmR_br$id = rownames(seurat_pheno_scran_norm_batchlmR_br)
rownames(seurat_pheno_scran_norm_batchlmR_br) <- exprs_batch_corrected_scran_colname_seurat
write.table(seurat_pheno_scran_norm_batchlmR_br, file=paste0(output_dir,prefix,".seurat.pheno.xls"),sep="\t",row.names=T,col.names=T,quote=F)
```


```{r scran-norm-batchlmR-seurat-initial, include=FALSE}
seurat_ct_scran_norm_batchlmR <- new('seurat', raw.data=exprs_batch_corrected_scran_seurat)
seurat_ct_scran_norm_batchlmR <- Seurat::Setup(seurat_ct_scran_norm_batchlmR, project=prefix, min.cells=1, min.genes=1, names.field=2, names.delim="@", is.expr=0, do.logNormalize=F, meta.data=seurat_pheno_scran_norm_batchlmR_br) 
```

识别在样品之间变化差异显著的基因用于样品聚类分析 (Figure \@ref(fig:scran-norm-batchlmR-mean-var-genes-seurat-plot))。首先计算每个基因的表达均值和离散度（标准差）；然后把基因按照表达值分为20个区间，计算每个区间内方差的Z-score，最小化低表达基因的检测波动引入的基因表达值的变化。

```{asis, ehco=FALSE}
Identifies genes that are outliers on a 'mean variability plot'. First, uses a function to calculate average expression (fxn.x) and dispersion (fxn.y) for each gene. Next, divides genes into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin. The purpose of this is to identify variable genes while controlling for the strong relationship between variability and average expression. 
```

```{r scran-norm-batchlmR-mean-var-genes-seurat-compute, include=FALSE, fig.cap="Mean-variability plot for all genes. Highly variable genes were labelled with their names."}
seurat_ct_scran_norm_batchlmR <- Seurat::MeanVarPlot(seurat_ct_scran_norm_batchlmR,y.cutoff = 1,x.low.cutoff = 1,fxn.x = Seurat::expMean,fxn.y = Seurat::logVarDivMean, x.high.cutoff = 15, do.spike = ercc)
#seurat_var_genes <- seurat_ct_scran_norm_batchlmR@var.genes
```

```{r scran-norm-batchlmR-mean-var-genes-seurat-plot, fig.cap="Mean-variability plot for all genes. Highly variable genes were labelled with their names. The X-axis represents the mean expression level, and for Y-axis represents the log(Variance/mean)."}
seurat_ct_scran_norm_batchlmR <- Seurat::MeanVarPlot(seurat_ct_scran_norm_batchlmR,y.cutoff = 1,x.low.cutoff = 1,fxn.x = Seurat::expMean,fxn.y = Seurat::logVarDivMean, x.high.cutoff = 15, do.spike = ercc, do.recalc = FALSE)
seurat_var_genes_scran_norm_batchlmR <- seurat_ct_scran_norm_batchlmR@var.genes
```

共鉴定出`r I(length(seurat_var_genes_scran_norm_batchlmR))`个样品间差异大的基因，用于区分细胞的基因表达异质性，选取其中四个绘制如 Figure \@ref(fig:scran-norm-batchlmR-mean-var-genes-seurat-4)。该步筛选表达差异大的基因时未考虑样品的分组信息，部分差异大的基因可能是“组内”差异大，也有可能组间差异大，对于样品的再分类会很有帮助。

```{r scran-norm-batchlmR-mean-var-genes-seurat-4, fig.cap="Violin-plot for 4 highly variable genes. Y-axis represents log-transformed scran (counts per million)."}
#print(paste0("There are ", length(seurat_var_genes_scran_norm_batchlmR), "highly variable genes identified."))
Seurat::VlnPlot(seurat_ct_scran_norm_batchlmR, seurat_var_genes_scran_norm_batchlmR[1:4], size.x.use = 5, size.y.use=9, size.title.use = 12, size.use = 0.8)
```

利用PCA方法根据高表达变化基因对样品进行聚类分组，结果如Figure \@ref(fig:scran-norm-batchlmR-pca-plot-based-on-seurat-var-genes)所示。

```{r scran-norm-batchlmR-pca-plot-based-on-seurat-var-genes, fig.height=7, fig.cap="PCA plot of sample correlation based on the most variable genes. Each point represents one sample. Different point colors represent different sample groups. Different point shape represents different sequencing batches."}
# Default using object@var.genes for PCA input
seurat_ct_scran_norm_batchlmR <- Seurat::PCA(seurat_ct_scran_norm_batchlmR, do.print=F)
pca_dim_first_scran_norm_batchlmR <- dim(seurat_ct_scran_norm_batchlmR@pca.x)[2]
Seurat::PCAPlot(seurat_ct_scran_norm_batchlmR,1,2,pt.size=2, pt.shape=pheno_data_col[2], do.return=TRUE) + scale_shape_manual(values=shapes)
```

Figure \@ref(fig:scran-norm-batchlmR-seurat-vizpca) 展示了与主成分1 (PC1)和主成分2 (PC2) 最相关的前40个基因。横轴代表基因在对应主成分的rank得分，纵轴为每个基因。在PC1中得分为正的基因，其表达值在Figure \@ref(fig:scran-norm-batchlmR-pca-plot-based-on-seurat-var-genes) 中右侧的样品中表达高；在PC1中得分为负的基因，其表达值在Figure \@ref(fig:scran-norm-batchlmR-pca-plot-based-on-seurat-var-genes) 中左侧的样品中表达高；在PC2中得分为正的基因，其表达值在Figure \@ref(fig:scran-norm-batchlmR-pca-plot-based-on-seurat-var-genes) 中上部分的样品中表达高；在PC2中得分为负的基因，其表达值在Figure \@ref(fig:scran-norm-batchlmR-pca-plot-based-on-seurat-var-genes) 中下部分的样品中表达高；


```{r scran-norm-batchlmR-seurat-vizpca, fig.height=7, fig.width=7, fig.cap="Visualize top 40 genes associated with principal components. These genes are thought as classification markers for samples. Genes positively contributed to principal component 1 (PC1) highly expressed in samples at the right part of PCA plot. Genes positively controbuted to principal component2 (PC2) highly expressed in samples at the top part of PCA plot."}
Seurat::VizPCA(seurat_ct_scran_norm_batchlmR,1:2, num.genes=40, font.size=0.7)
```

Figure \@ref(fig:scran-norm-batchlmR-heatmap-pc1-seurat) 展示了6个主成分的关键基因在样品中的分布。每一列代表一个样品，横轴代表一个基因。样品和基因都根据其主成分得分排序，不同的子图的样品顺序不一定一致。颜色从红色-到黑色-到黄色表示基因的相对表达量由低到高。

```{r, eval=FALSE}
a = Seurat::PCHeatmap(seurat_ct_scran_norm_batchlmR, pc.use=1:6, do.balanced=T, label.columns = FALSE)
```

```{r scran-norm-batchlmR-heatmap-pc1-seurat, fig.height=7, fig.cap="Heatmap showing sources of heterogeneity for PC1-6."}
Seurat::PCHeatmap(seurat_ct_scran_norm_batchlmR,pc.use=1:6, do.balanced=T, label.columns = FALSE, remove.key = FALSE)
```
```{r scran-Determine-statistically-significant-principal-components, eval=F}
# Do 200 random samplings to find significant genes, each time randomly permute 1% of genes
# This returns a 'p-value' for each gene in each PC, based on how likely the gene/PC score woud have been observed by chance
# Note that in this case we get the same result with 200 or 1000 samplings, so we do 200 here for expediency
seurat_ct_scran_norm_batchlmR=Seurat::JackStraw(seurat_ct_scran_norm_batchlmR,num.replicate = 200,do.print = FALSE)

# The jackStraw plot compares the distribution of P-values for each PC with a uniform distribution (dashed line)
# 'Significant' PCs will have a strong enrichment of genes with low p-values (solid curve above dashed line)
# In this case PC1-9 are strongly significant
Seurat::JackStrawPlot(seurat_ct_scran_norm_batchlmR,PCs = 1:12)
```

```{r scran-norm-batchlmR-grow-gene-list-pca-seurat, fig.height=7, eval=FALSE}
# Previous analysis was performed on < 400 variable genes. To identify a larger gene set that may drive
# biological differences, but did not pass our mean/variability thresholds, we first calculate PCA scores
# for all genes (PCA projection)
seurat_ct_scran_norm_batchlmR <- Seurat::ProjectPCA(seurat_ct_scran_norm_batchlmR, do.print=F, do.center=T)
# For test
# Seurat::PrintPCA(seurat_ct_scran_norm_batchlmR,1)
# Seurat::PCHeatmap(seurat_ct_scran_norm_batchlmR,pc.use=1, use.full=TRUE, do.balanced=T, remove.key = T, num.genes=40)
```

```{r scran-norm-batchlmR-old-seirat-sig-genes, eval=F}
# Choose significant genes for PC1-9, allow each PC to contribute a max of 200 genes (to avoid one PC swamping the analysis)
seurat_ct_scran_norm_batchlmR_sig_genes = Seurat::PCASigGenes(seurat_ct_scran_norm_batchlmR, pcs.use=1:pca_dim_first_scran_norm_batchlmR, pval.cut=0.01, max.per.pc = 200)
length(seurat_ct_scran_norm_batchlmR_sig_genes)
```

```{r, eval=F}
# Now redo the PCA analysis with the new gene list
seurat_ct_scran_norm_batchlmR <- Seurat::PCA(seurat_ct_scran_norm_batchlmR, pc.genes=seurat_ct_scran_norm_batchlmR_sig_genes, do.print=FALSE)

# Redo random sampling, PCs 1-11 are significant (additional PCs are borderline, but signal disappears with 1k replicates, though we do 200 here for expediency)
seurat_ct_scran_norm_batchlmR =Seurat::JackStraw(seurat_ct_scran_norm_batchlmR,num.replicate = 200,do.print = FALSE)
Seurat::JackStrawPlot(seurat_ct_scran_norm_batchlmR, PCs = 1:15)
```

```{asis, echo=FALSE}
PC selection – identifying the true dimensionality of a dataset – is an important step for Seurat, but can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-10 as a cutoff. We followed the jackStraw here, admittedly buoyed by seeing the PCHeatmap returning interpretable signals (including canonical dendritic cell markers) throughout these PCs. Though the results are only subtly affected by small shifts in this cutoff (you can test below), we strongly suggest always explore the PCs they choose to include downstream.
```

```{r scran-norm-batchlmR-seurat-tsne-plot, fig.height=10, fig.width=8, fig.cap="tSNE map of normalized expression data with batch effect removed given different perplexity values using Seurat. Multiple iterations (1000) were performed to make sure the stability of clustering.", eval=F}
pca_dim_scran_norm_batchlmR <- dim(seurat_ct_scran_norm_batchlmR@pca.x)[2]
plot_tsne_list_scran_norm_batchlmR <- by(data=cluster_perplexity, INDICES = cluster_perplexity$a, FUN= function(x) {
  x <- x$a
  seurat_ct_scran_norm_batchlmR <- Seurat::RunTSNE(seurat_ct_scran_norm_batchlmR, dims.use=1:pca_dim_scran_norm_batchlmR, k.seed=11521, max_iter=1000, perplexity = x)
  Seurat::TSNEPlot(seurat_ct_scran_norm_batchlmR, pt.size=1, do.return=T, no.legend=TRUE, do.label=T) +
    labs(title=paste0("perplexity = ",x)) + theme(legend.text=element_text(size=5))
})

do.call(grid.arrange, c(plot_tsne_list_scran_norm_batchlmR, ncol=2))
```

```{r scran-norm-batchlmR-seurat-regroup-tens-test, fig.width=7, fig.height=7, fig.cap="Test different parameters for density clustering to make sure get the best re-classification.", eval=F}
perplexity= 4
#seurat_ct_scran_norm_batchlmR <- Seurat::RunTSNE(seurat_ct_scran_norm_batchlmR, dims.use=1:pca_dim_scran_norm_batchlmR, k.seed=11521, max_iter=2000, perplexity = perplexity)

plot_tsne_list <- by(data=cluster_perplexity, INDICES = cluster_perplexity$a, FUN= function(x) {
  x <- x$a
  seurat_ct_scran_norm_batchlmR <- Seurat::DBClustDimension(seurat_ct_scran_norm_batchlmR,dim.1=1,dim.2=2,reduction.use = "tsne", G.use=x, set.ident=T)
  a = Seurat::TSNEPlot(seurat_ct_scran_norm_batchlmR, pt.size=1, pt.shape=pheno_data_col[2], do.return=T, no.legend=TRUE, do.label=T) +
    labs(title=paste0("G.use = ",x)) + theme(legend.text=element_text(size=5))
  #seurat_ct_scran_norm_batchlmR = Seurat::SetAllIdent(seurat_ct_scran_norm_batchlmR, "orig.ident")
  a
})

do.call(grid.arrange, c(plot_tsne_list, ncol=2))
```

```{r scran-seurat-regroup-tens, fig.width=7, fig.height=7, fig.cap="", eval=F}
#G_use = 4 
#seurat_ct_scran_norm_batchlmR <- Seurat::DBClustDimension(seurat_ct_scran_norm_batchlmR,dim.1=1,dim.2=2,reduction.use = "tsne", G.use=G_use, set.ident=T)
#Seurat::TSNEPlot(seurat_ct_scran_norm_batchlmR, pt.size=1, pt.shape=pheno_data_col[1], do.return=T) +
#  labs(title=paste0("G.use = ",G_use)) + theme(legend.text=element_text(size=5))

#seurat_identity = Seurat::FetchData(seurat_ct_scran_norm_batchlmR, c("ident", "orig.ident"))
#write.table(seurat_identity, file=paste0(dir,prefix,".seurat_recalssification_phenoData.xls"), sep="\t", row.names=T, col.names=T, quote=F)
## You can also switch easily switch the cell's identity (for example, going back to the original annotation)
#nbt = set.all.ident(nbt, "orig.ident")
## And switch back - to the cluster ID defined by the tree building
#nbt = set.all.ident(nbt, "tree.ident")
```

```{r scran-buildClusterTree-seurat, eval=F}
# Build a phylogenetic tree, and rename/reorder cluster names according to their position on the tree
# See help for details on tree building strategy
# This gives closely related clusters similar cluster IDs, which is occasionally useful for visualization later on
# Assigned cluster will be placed in the 'tree.ident' field of nbt@data.info, and also stored in nbt@ident
#seurat_ct_scran_norm_batchlmR = Seurat::BuildClusterTree(seurat_ct_scran_norm_batchlmR, do.reorder = TRUE, reorder.numeric = TRUE,pcs.use = 1:pca_dim_scran_norm_batchlmR)
#seurat_identity = Seurat::FetchData(seurat_ct_scran_norm_batchlmR, c("ident", "orig.ident"))
#write.table(seurat_identity, file=paste0(dir,prefix,".seurat_recalssification_phenoData.xls"), sep="\t", #row.names=T, col.names=T, quote=F)
## You can also switch easily switch the cell's identity (for example, going back to the original annotation)
#nbt = set.all.ident(nbt, "orig.ident")
## And switch back - to the cluster ID defined by the tree building
#nbt = set.all.ident(nbt, "tree.ident")
```

```{r scran-buildClusterTree-seurat-label-tsne, fig.cap="", eval=FALSE}
#Seurat::TSNEPlot(seurat_ct_scran_norm_batchlmR, pt.size=1, pt.shape=pheno_data_col[1], do.return=T, do.label=T, label.pt.size=0.5) +
#  labs(title=paste0("Re-classsification")) + theme(legend.text=element_text(size=5))
```

Figure \@ref(fig:scran-norm-batchlmR-stanrarddeviation-pc) 展示了不同主成分对样品差异的贡献，便于筛选主成分用于后续分析。横轴代表识别出的主成分，1、2、3...。纵轴代表主成分反应的样品差异的度量，数值越大，对应的主成分反映的样品差异越多。

```{r scran-norm-batchlmR-stanrarddeviation-pc, fig.cap="Variations contributed by each principal component. "}
Seurat::PCElbowPlot(seurat_ct_scran_norm_batchlmR)
```

Figure \@ref(fig:scran-norm-batchlmR-seurat-findclusters) and \@ref(fig:scran-norm-batchlmR-seurat-findclusters-0304) 展示了根据样品表达图谱进行的不同精细程度的分类。resolution越大，分出的类越多。可以直接选择大的resolution分成小类，也可以先分出大类，然后再细分亚类。

```{r scran-norm-batchlmR-seurat-findclusters, fig.height=13, fig.cap="Test different parameters for density clustering to make sure get the best re-classification."}
findc_resolution_scran_norm_batchlmR <- data.frame(res=c(0.4,0.7,1.0,1.5,2,4))
pca_dim_scran_norm_batchlmR <- dim(seurat_ct_scran_norm_batchlmR@pca.x)[2]

find_cluster_tsne_list_scran_norm_batchlmR <- by(data=findc_resolution_scran_norm_batchlmR, INDICES = findc_resolution_scran_norm_batchlmR$res, FUN=function(x){
  res <- x$res
  seurat_ct_scran_norm_batchlmR <- Seurat::FindClusters(seurat_ct_scran_norm_batchlmR, pc.use=1:pca_dim_scran_norm_batchlmR, resolution=res, print.output=0, random.seed=11521, save.SNN=T, temp.file.location="./")
  seurat_ct_scran_norm_batchlmR <- Seurat::RunTSNE(seurat_ct_scran_norm_batchlmR, dims.use=1:pca_dim_scran_norm_batchlmR, k.seed=11521, do.fast=T, perplexity = round(sampleNumber_scran/5))
  Seurat::TSNEPlot(seurat_ct_scran_norm_batchlmR, pt.size=1, do.return=T, do.label=T, no.legend=T) +
    labs(title=paste0("perplexity = ",round(sampleNumber_scran/5), "; resolution = ", res)) +
    theme(legend.text=element_text(size=5))
})


do.call(grid.arrange, c(find_cluster_tsne_list_scran_norm_batchlmR, ncol=2))
```

```{r scran-norm-batchlmR-seurat-findclusters-0304, fig.height=13, fig.cap="Test different parameters for density clustering to make sure get the best re-classification."}

find_cluster_tsne_list_scran_norm_batchlmR0304 <- by(data=findc_resolution_scran_norm_batchlmR, INDICES = findc_resolution_scran_norm_batchlmR$res, FUN=function(x){
  res <- x$res
  seurat_ct_scran_norm_batchlmR <- Seurat::FindClusters(seurat_ct_scran_norm_batchlmR, pc.use=1:pca_dim_scran_norm_batchlmR, resolution=res, print.output=0, random.seed=304, save.SNN=T, temp.file.location="./")
  seurat_ct_scran_norm_batchlmR <- Seurat::RunTSNE(seurat_ct_scran_norm_batchlmR, dims.use=1:pca_dim_scran_norm_batchlmR, k.seed=304, do.fast=T, perplexity = round(sampleNumber_scran/5))
  Seurat::TSNEPlot(seurat_ct_scran_norm_batchlmR, pt.size=1, do.return=T, do.label=T, no.legend=T) +
    labs(title=paste0("perplexity = ",round(sampleNumber_scran/5), "; resolution = ", res)) +
    theme(legend.text=element_text(size=5))
})


do.call(grid.arrange, c(find_cluster_tsne_list_scran_norm_batchlmR0304, ncol=2))
```



```{r scran-norm-batchlmR-seurat-findclusters-final, fig.cap="Sample re-classification."}
res_scran_norm_batchlmR=1.6
seed_scran_norm_batchlmR=11521
seurat_ct_scran_norm_batchlmR <- Seurat::FindClusters(seurat_ct_scran_norm_batchlmR, pc.use=1:pca_dim_scran_norm_batchlmR, resolution=res_scran_norm_batchlmR, print.output=0, random.seed=seed_scran_norm_batchlmR, save.SNN=T, temp.file.location="./")
seurat_ct_scran_norm_batchlmR <- Seurat::RunTSNE(seurat_ct_scran_norm_batchlmR, dims.use=1:pca_dim_scran_norm_batchlmR, k.seed=seed_scran_norm_batchlmR, do.fast=T, perplexity = round(sampleNumber_scran/5))
Seurat::TSNEPlot(seurat_ct_scran_norm_batchlmR, pt.size=1, do.return=T, do.label=T, no.legend=TRUE) +
    labs(title=paste0("perplexity = ",round(sampleNumber_scran/5), "; resolution = ", res_scran_norm_batchlmR)) +
    theme(legend.text=element_text(size=5))
seurat_identity_scran_norm_batchlmR = Seurat::FetchData(seurat_ct_scran_norm_batchlmR, c("ident", "orig.ident"))
write.table(seurat_identity_scran_norm_batchlmR, file=paste0(output_dir,prefix,"seurat_scran_norm_log2_batchlmR_batch_rmove_recalssification_phenoData.xls"), sep="\t", row.names=T, col.names=T, quote=F)
```


选取**resolution=`r res_scran_norm_batchlmR`**作为分类标准，进行下游探索 (Figure \@ref(fig:scran-norm-batchlmR-seurat-findclusters-final))。分类结果以及新分类与原分类的对应表格可[点击下载](`r  I(paste0("scater/", prefix,".seurat_scran_norm_batchlmR_lm_batch_rmove_recalssification_phenoData.xls"))`)。


```{r scran-norm-batchlmR-save-seurat-reclassification}
save(seurat_ct_scran_norm_batchlmR,file=paste0(prefix, ".seurat_scran_norm_log2_batchlmR_batch_rmove_reclassification.Robj"))
```


```{r scran-norm-batchlmR-find-all-markers}
# thresh.use: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. 
# min.pct: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. 
all_seurat_markers_scran_norm_batchlmR = Seurat::FindAllMarkers(seurat_ct_scran_norm_batchlmR, thresh.use=0.25, min.pct=0.2, test.use="roc", print.bar = debug, only.pos=TRUE)
ident_scran_norm_batchlmR <- unique(seurat_identity_scran_norm_batchlmR$ident)
```


样品标记基因列表 (Table \@ref(tab:scran-norm-batchlmR-list-marker-genes))。

Table: (\#tab:scran-norm-batchlmR-marker-gene-label-explanation) Explanation for each column of marker files.

```{r scran-norm-batchlmR-marker-gene-label-explanation, results="asis"}
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

**样品标记基因示例**

```{r scran-norm-batchlmR-sample-marker-gene-example-top5}
sample_marker_gene_exp_top5_scran_norm_batchlmR <- all_seurat_markers_scran_norm_batchlmR %>% group_by(cluster) %>% top_n(5, avg_diff)
knitr::kable(sample_marker_gene_exp_top5_scran_norm_batchlmR, booktabs=T, caption="List of top 5 markers for each cluster.")
```

Table: (\#tab:scran-norm-batchlmR-list-marker-genes) List of markers for each sample class. Currently only positive genes are reported to shrink computation time.

```{r scran-norm-batchlmR-list-marker-genes, results="asis"}
ident2_scran_norm_batchlmR <- ident_scran_norm_batchlmR[order(ident_scran_norm_batchlmR)]

listMarkers_scran_norm_batchlmR <- function(cluster){
  c(paste0("[Cluster ", cluster," UP marker](",paste0("scater/", prefix,".seurat.scran_norm_log2transform_batchlmR_batch_rmove.",cluster,'.up_marker.xls'),")"),
  paste0("[Cluster ", cluster," DW marker](",paste0("scater/", prefix,".seurat.scran_norm_log2transform_batchlmR_batch_rmove.",cluster,'.dw_marker.xls'),")"))
}
marker_list_tmp_scran_norm_batchlmR <- as.data.frame(t(sapply(ident2_scran_norm_batchlmR, listMarkers_scran_norm_batchlmR)))
colnames(marker_list_tmp_scran_norm_batchlmR) <- c("UP","DW (ignore)")
rownames(marker_list_tmp_scran_norm_batchlmR) <- ident2_scran_norm_batchlmR
knitr::kable(marker_list_tmp_scran_norm_batchlmR, booktabs=T, format="markdown", caption="List of markers for each sample class.")
```


```{r scran-norm-batchlmR-output-all-markers}
outputMarkers_scran_norm_batchlmR <- function(cluster){
  up_marker <- all_seurat_markers_scran_norm_batchlmR[all_seurat_markers_scran_norm_batchlmR$cluster==cluster & all_seurat_markers_scran_norm_batchlmR$avg_diff>0,]
  write.table(up_marker, file=paste0(output_dir,prefix,".seurat.scran_norm_log2transform_batchlmR_batch_rmove.",cluster,'.up_marker.xls'), sep="\t", row.names=F, col.names=T, quote=F)
  dw_marker <- all_seurat_markers_scran_norm_batchlmR[all_seurat_markers_scran_norm_batchlmR$cluster==cluster & all_seurat_markers_scran_norm_batchlmR$avg_diff<0,]
  write.table(dw_marker, file=paste0(output_dir,prefix,".seurat.scran_norm_log2transform_batchlmR_batch_rmove.",cluster,'.dw_marker.xls'), sep="\t", row.names=F, col.names=T, quote=F)
  c(up_marker$gene[1], up_marker$gene[2])
}
test_marker_for_plot_scran_norm_batchlmR <- sapply(ident_scran_norm_batchlmR, outputMarkers_scran_norm_batchlmR)
```

Figure \@ref(fig:scran-norm-batchlmR-plot-representative-markers) 展示了代表性样品Marker基因在对应样品的表达图谱。

```{r scran-norm-batchlmR-plot-representative-markers, fig.height=18, fig.cap="Violin-plot for representative markers. Y-axis represents log-transformed scran (counts per million)."}
test_marker_for_plot_scran_norm_batchlmR <- test_marker_for_plot_scran_norm_batchlmR[!is.na(test_marker_for_plot_scran_norm_batchlmR)]
Seurat::VlnPlot(seurat_ct_scran_norm_batchlmR, test_marker_for_plot_scran_norm_batchlmR, size.x.use = 5, size.y.use=6, size.title.use = 10, size.use = 0.8, nCol=2)
```

Figure \@ref(fig:scran-norm-batchlmR-seurat-featureplot) 展示数个样品标记基因在样品间的表达情况（灰色到蓝色表达由低到高）。

```{r scran-norm-batchlmR-seurat-featureplot, fig.height=7, fig.cap="Map of marker genes to related sample classifications. At most two marker genes were selected for each sample.  Different colors represent the expression (low:grey;high:blue) of the gene in the context of all the cells and help validate the specificity of the marker."}
Seurat::FeaturePlot(seurat_ct_scran_norm_batchlmR, test_marker_for_plot_scran_norm_batchlmR,cols.use = c("grey","dark blue"))
```

Figure \@ref(fig:scran-norm-batchlmR-seurat-heatmap-20-markers) 展示了每个cluster最具代表性的10个基因的相对表达图谱 （颜色从红色-黑色-黄色代表相对表达值由低到高）。

```{r scran-norm-batchlmR-seurat-heatmap-20-markers, fig.height=14, fig.cap="Heatmap showing expression patterns of top 20 marker genes for each group. Colors from read to black to yellow showing Z-scaled expression value from low to high."}
all_seurat_markers_scran_norm_batchlmR %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_diff) -> top10_scran_norm_batchlmR
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
Seurat::DoHeatmap(seurat_ct_scran_norm_batchlmR, genes.use = top10_scran_norm_batchlmR$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)
```

END
