

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Seurat)
library(tidyverse)
library(NMF)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(WGCNA)
library(reshape2)
library(stringr)
library(SingleR)

setwd("D:/shangke/lession17/")
load("scRNA_harmony.Rdata")
load("D:/ref.data/ref_Human_all.RData")

refdata <- ref_Human_all

testdata <- GetAssayData(scRNA_harmony, slot="data")
###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
clusters <- scRNA_harmony@meta.data$seurat_clusters
###开始用singler分析
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA_harmony@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(scRNA_harmony@meta.data$celltype)
table(scRNA_harmony@meta.data$seurat_clusters)


scRNA_harmony=scRNA_harmony[,sample(colnames(scRNA_harmony),1500)] 




setwd("D:/shangke/lession25/")
 
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData(do.center = F)

## 使用pca的分解结果降维聚类
scRNA_harmony <- RunPCA(scRNA_harmony)
set.seed(803)
scRNA_harmony.pca <- RunUMAP(scRNA_harmony, dims = 1:10) %>% FindNeighbors(dims = 1:10) %>% FindClusters()

## 结果可视化
  DimPlot(scRNA_harmony.pca, label = T) + ggsci::scale_color_igv()
 
 

############################################################################################################

## 高变基因表达矩阵的分解
#  
vm <- scRNA_harmony@assays$RNA@scale.data


NMF包通过nmf()函数实现矩阵分解，它的用法及重要参数如下：
nmf(x, rank, method, seed, nrun, ...)

x：待分解非负矩阵，数据格式可以是matrix，data.frame， ExpressionSet
rank：分解的基数量，对于单细胞数据，可以设置为期望的细胞类型数量或表达模式数量
method：因式分解的常用方法，这里介绍三种常用的
1、基于KL 散度进行度量目标函数的多重迭代梯度下降算法——brunet(默认算法)
2、基于欧几里得距离度量目标函数的多重迭代梯度下降算法——lee
3、交替最小二乘法(Alternating Least Squares(ALS))——snmf/r  
seed：因式分解的初始化种子    
nrun：运行次数       


res <- nmf(vm, 10, method = "snmf/r", seed = 'nndsvd') 
runtime(res)
 

## 分解结果返回suerat对象
scRNA_harmony@reductions$nmf <- scRNA_harmony@reductions$pca
scRNA_harmony@reductions$nmf@cell.embeddings <- t(coef(res))    
scRNA_harmony@reductions$nmf@feature.loadings <- basis(res)  

## 使用nmf的分解结果降维聚类
set.seed(803)
scRNA_harmony.nmf <- RunUMAP(scRNA_harmony, reduction = 'nmf', dims = 1:10) %>% 
  FindNeighbors(reduction = 'nmf', dims = 1:10) %>% FindClusters()

## 结果可视化  
 DimPlot(scRNA_harmony.nmf, label = T) + ggsci::scale_color_igv()
 
 
 ## 查看细胞因子上的荷载
 tmp <- data.frame(t(coef(res)), check.names = F)
 colnames(tmp) <- paste0("factor", 1:10)
 scRNA_harmony.nmf <- AddMetaData(scRNA_harmony.nmf, metadata = tmp)
  FeaturePlot(scRNA_harmony.nmf, features = paste0("factor", 1:9), ncol = 3)
 
 
 ## 查看细胞主成分上的荷载
  FeaturePlot(scRNA_harmony.nmf, features = paste0("PC_", 1:9), ncol = 3)

  对比上下两张图，很容易发现NMF的因子比PCA的PC轴解释性更强。  

  "COL1A1","CDH11","RUNX2","CTSK","MMP9" ,"PCNA"
  
  FeaturePlot(scRNA_harmony.nmf, features = c( "COL1A1","CDH11","RUNX2","CTSK","MMP9" ,"PCNA"), ncol = 2)
  FeaturePlot(scRNA_harmony.pca, features = c( "COL1A1","CDH11","RUNX2","CTSK","MMP9" ,"PCNA"), ncol = 2)
  

  提取celltype的signatures
  
  ## 提取每个因子贡献度最大的20个基因
  f <- extractFeatures(res, 20L)
  f <- lapply(f, function(x) rownames(res)[x])
  f <- do.call("rbind", f)
  DT::datatable(t(f))

