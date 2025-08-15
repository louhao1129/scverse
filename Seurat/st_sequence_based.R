# 本教程演示了如何使用 Seurat（>=3.2）分析空间分辨的 RNA-seq 数据。虽然分析流程与单细胞 RNA-seq 分析的 Seurat 工作流程相似，
# 但我们介绍了更新的交互和可视化工具，特别强调空间和分子信息的整合。
# 本教程将涵盖以下任务，我们认为这些任务将对许多空间分析很常见：

# Normalization  归一化
# Dimensional reduction and clustering 降维和聚类
# Detecting spatially-variable features 检测空间变异性特征
# Interactive visualization 交互式可视化
# Integration with single-cell RNA-seq data 与单细胞 RNA 测序数据的整合
# Working with multiple slices 处理多个切片


# 我们分析了一个使用 10x Genomics 的 Visium 技术生成的数据集

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
# 在 Seurat 对象中，spot与基因表达矩阵类似于典型的“RNA” Assay ，但包含的是spot水平数据，而不是单细胞水平数据。
# 图像本身存储在 Seurat 对象的新的 images 槽中。 images 槽还存储了将斑点与其在组织图像上的物理位置相关联所需的信息。

# 数据预处理
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
