rm(list = ls())
gc()
setwd("./scverse/Signac")

# 文件准备1：
# 1.atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
# 2.atac_v1_pbmc_10k_singlecell.csv
# 3.atac_v1_pbmc_10k_fragments.tsv.gz

# 加载需要包----
# options(download.file.method = "curl")
# install.packages("Signac")
# install.packages("Seurat")
# BiocManager::install("biovizBase")
# BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86')) # 人 hg38
# BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'EnsDb.Hsapiens.v75')) # 人 hg19 
# BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79')) # 鼠mm10
# install.packages("BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz", repos = NULL, type = "source")
# install.packages("BSgenome.Mmusculus.UCSC.mm10_1.4.3.tar.gz", repos = NULL, type = "source")

library(Signac)
library(Seurat)
# library(EnsDb.Hsapiens.v86) # 人
library(EnsDb.Hsapiens.v75) # 人
# library(EnsDb.Mmusculus.v79) # 鼠
# library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
# library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
library(future)

# 分布式计算 适用与scRNA分析和scATAC
plan(multicore, workers = 3) # windows:multisession  linux:multisession,multicore
# system("wmic cpu get NumberOfCores, NumberOfLogicalProcessors") # windows查看系统核心数
options(future.globals.maxSize = 10 * 1024^3) # for 10 Gb RAM


path_10X_samples <- "./data/GSE000"
samples <- list.files(path_10X_samples)
samples <- samples[gtools::mixedorder(samples)]
atacS4_list <- list()
samples <- samples[2:3] # 选取两个样本
for(sample in samples){
  # sample = samples[1]
  data.path <- file.path(path_10X_samples, sample)
  files <- dir(data.path) # List the Files in a Directory/Folder
  peak_bc_matrix <- files[stringr::str_detect(files,"peak_bc_matrix.h5")]
  singlecell <- files[stringr::str_detect(files,"singlecell.csv")]
  frag <- files[stringr::str_detect(files,"fragments.tsv.gz$")]
  
  counts <- Read10X_h5(filename = file.path(data.path,peak_bc_matrix))
  metadata <- read.csv(
    file = file.path(data.path,singlecell),
    header = TRUE,
    row.names = 1
  )
  chrom_assay <- CreateChromatinAssay(
    counts = counts,  # 输入counts
    sep = c(":", "-"),
    fragments = file.path(data.path,frag),
    min.cells = 10,
    min.features = 200
  )
  scATAC_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  atacS4_list <- append(atacS4_list, scATAC_obj)
}

scATAC_object <- merge(atacS4_list[[1]], 
                       y = atacS4_list[-1],
                       add.cell.ids = samples)
# save(scATAC_object,file = "F:/sxdw/scATAC/class7/rdata/scATAC_object.rda")
# load("F:/sxdw/scATAC/class7/rdata/scATAC_object.rda")
scATAC_object$sample <- stringr::str_split(colnames(scATAC_object),"_")[[1]][1]

granges(scATAC_object)



# scATAC的S4对象添加ann和seqinfo信息 ---- 
# 提取基因信息
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # 人hg38
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75) # 人hg19
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) # 鼠
# 修改染色体名称 # 修改seqnames中的values
seqlevels(annotations) <- paste0('chr', seqlevels(annotations)) 

# 修改基因组名称 # 修改seqinfo中的genome 
genome(annotations) <- "hg38"  # 人
# genome(annotations) <- "hg19"  # 人
# genome(annotations) <- "mm10" # 鼠
# 添加注释信息到对象
Annotation(scATAC_object) <- annotations # 添加到assays - peaks - annotation
rm(list=ls()[which(!(ls() %in% c("scATAC_object")))])

# 计算QC指标 meta.data---- 
# 核小体信号评分：反应染色质开放状态 
scATAC_object <- NucleosomeSignal(object = scATAC_object)
scATAC_object$nucleosome_group <- ifelse(scATAC_object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# 添加黑名单比例信息和峰/reads比例
# scATAC_object$pct_reads_in_DNase <- scATAC_object$DNase_sensitive_region_fragments / scATAC_object$passed_filters * 100
scATAC_object$pct_reads_in_peaks <- scATAC_object$peak_region_fragments / scATAC_object$passed_filters * 100
scATAC_object$blacklist_ratio <- scATAC_object$blacklist_region_fragments / scATAC_object$peak_region_fragments

# 转录起始位点富集评分：反应细胞转录活动信号 # peaks and meta.data
# 当fast=T时，仅计算TSS分数，而不存储每个细胞Tn5插入频率的位置矩阵，这样后续将不能使用TSSPlot画图。
scATAC_object <- TSSEnrichment(object = scATAC_object, fast = FALSE)

# 用于快速找到不同QC的合适截止值
# pdf("计算QC指标.pdf",width = 6, height = 6)
DensityScatter(scATAC_object, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
# dev.off()

scATAC_object$high.tss <- ifelse(scATAC_object$TSS.enrichment > 3, 'High', 'Low')




pdf("nucleosome_group.pdf",width = 6, height = 5)
FragmentHistogram(object = scATAC_object, group.by = 'nucleosome_group')
dev.off()

# 基于TSS分数将不同的细胞进行分组以及画所有TSS位点的可及性信号图来检查TSS富集分数
pdf("high.tss.pdf",width = 6, height = 5)
TSSPlot(scATAC_object, group.by = 'high.tss') + NoLegend()
dev.off()

pdf("指标汇总.pdf",width = 8, height = 8)
VlnPlot(
  object = scATAC_object,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 3
)
dev.off()

# QC ---- 
dim(scATAC_object)
scATAC_object <- subset(
  x = scATAC_object,
  subset = nCount_peaks > 1000 & # 1000
    nCount_peaks < 30000 & # 50000
    # pct_reads_in_DNase > 40 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
dim(scATAC_object)


# 归一化和降维----

# LSI分析:TF-IDF 和 SVD 的组合步骤z
# TF-IDF归一化。这是一个两步归一化过程，既可以在细胞间进行归一化以校正细胞测序深度的差异，又可以在峰之间进行归一化，以便为更罕见的峰提供更高的值。
scATAC_object <- RunTFIDF(scATAC_object)
# 寻找高变峰 “q75”表示使用前25%的峰 储存在 VariableFeatures
scATAC_object <- FindTopFeatures(scATAC_object, min.cutoff = 'q0')
# SVD矩阵分解，类似于PCA的输出
scATAC_object <- RunSVD(scATAC_object)

pdf("每个 LSI 分量与测序深度之间的相关性.pdf",height = 4,width = 8)
DepthCor(scATAC_object,n=NULL)
dev.off()


scATAC_object <- RunUMAP(object = scATAC_object, reduction = 'lsi', dims = 2:50)
scATAC_object <- FindNeighbors(object = scATAC_object, reduction = 'lsi', dims = 2:50)
scATAC_object <- FindClusters(object = scATAC_object, verbose = FALSE, algorithm = 3,resolution = 0.8)

pdf("umap.pdf",height = 5,width = 5)
DimPlot(object = scATAC_object, label = TRUE) + NoLegend()
dev.off()

save(scATAC_object,file = "F:/sxdw/scATAC/class7/rdata/class7.rda")

