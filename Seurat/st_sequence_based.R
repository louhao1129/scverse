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

getwd()
setwd("./scverse/Seurat")

# InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
# 在 Seurat 对象中，spot与基因表达矩阵类似于典型的“RNA” Assay ，但包含的是spot水平数据，而不是单细胞水平数据。
# 图像本身存储在 Seurat 对象的新的 images 槽中。 images 槽还存储了将spot与其在组织图像上的物理位置相关联所需的信息。

# 数据预处理
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# 这些图表明，spot之间分子计数的方差不仅具有技术性，还取决于组织解剖结构。例如，组织中对神经元缺乏的区域（如皮质白质），可重复地表现出较低的分子计数。
# 因此，标准方法（如 LogNormalize() 函数），在标准化后强制使每个数据点具有相同的潜在“大小”，可能会出现问题。
# 作为替代方案，我们推荐使用 sctransform

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# 基因表达可视化
# SpatialFeaturePlot() 函数扩展了 FeaturePlot() ，可以在组织病理学上叠加分子数据
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

plot <- SpatialFeaturePlot(brain, features = c("Ttr")) + theme(legend.text = element_text(size = 0),
    legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
jpeg(filename = "./output/images/spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()
plot

# 降维、聚类和可视化
# 我们可以接着对 RNA 表达数据运行降维和聚类，使用与 scRNA-seq 分析相同的流程。
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# 我们可以在 UMAP 空间（使用 DimPlot() ）中可视化聚类结果，或将其叠加在图像上（使用 SpatialDimPlot() ）。
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

# 交互式绘图
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain)

# 空间可变特征的识别

# Seurat 提供两种工作流程来识别与组织中空间位置相关的分子特征。第一种是基于组织中预先注释的解剖区域进行差异表达分析，这些区域可能通过无监督聚类或先验知识确定。
# 在这个案例，这种策略效果很好，因为上述聚类表现出明显的空间限制。

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

# 种替代方法，在 FindSpatiallyVariables() 中实现，是在没有预先注释的情况下搜索表现出空间模式性的特征。该过程计算 gamma(r) 值，用于测量两个相距特定“r”距离的点之间的依赖关系
# 献中有多种方法可以完成这项任务，包括 SpatialDE 和 Splotch。
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],selection.method = "moransi")

# 可视化由该度量识别的前 6 个特征的表达情况
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

# 提取解剖区域
# 与单细胞对象类似，您可以对对象进行子集化以专注于数据子集


# 与单细胞数据的整合
# 此处介绍的标签转移工作流，利用 sctransform 标准化