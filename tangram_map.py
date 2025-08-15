import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import os
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg
# sc.logging.print_header()
print(f"squidpy=={sq.__version__}")
print(os.getcwd())
os.chdir("./scverse/")
# 加载数据集
# 加载 Squidpy 中可用的公共数据，来源为小鼠脑皮层。
# 单细胞数据存储在 adata_sc 中，空间数据存储在 adata_st 中
# adata_st = sq.datasets.visium_fluo_adata_crop()
adata_st = sc.read_h5ad("./data/tangram_tutorial/tangram_st.h5ad")
adata_st = adata_st[
    adata_st.obs.cluster.isin([f"Cortex_{i}" for i in np.arange(1, 5)])
].copy()
# img = sq.datasets.visium_fluo_image_crop()
# 我们将小鼠脑的裁剪部分仅包含脑皮层的簇
# adata_sc = sq.datasets.sc_mouse_cortex()
adata_sc = sc.read_h5ad("./data/tangram_tutorial/tangram_sc.h5ad")

# 可视化空间和单细胞数据集
fig, axs = plt.subplots(1, 2, figsize=(20, 5))
sc.pl.spatial(
    adata_st, color="cluster", alpha=0.7, frameon=False, show=False, ax=axs[0]
)
sc.pl.umap(
    adata_sc, color="cell_subclass", size=10, frameon=False, show=False, ax=axs[1]
)
# plt.tight_layout()

# Tangram 通过观察用户指定的训练基因（即一组基因）来学习单个细胞数据的空间对齐。训练基因需要具有有趣的信号并具有高质量的测量值。
# 通常，我们选择 100-1000 个差异表达基因作为训练基因，这些基因按细胞类型分层。
# 有时，我们也使用整个转录组，或使用不同的训练基因集进行不同的映射，以查看结果的变化程度。
# Tangram 使用基于余弦相似度的自定义损失函数，将 scRNA-seq 谱图拟合到空间中

# 预处理
# 在这里例子里面，我们使用 1401 个标记基因作为训练基因。(标记基因选取为每群细胞特异
# 表达的前100个基因）
sc.tl.rank_genes_groups(adata_sc, groupby="cell_subclass", use_raw=False) # use_raw=False 指定不使用原始（未处理的）数据进行分析。
markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
len(markers)

# 我们使用 pp_adatas 准备数据，它执行以下操作：
# 通过 genes 参数从用户处获取一组基因。这些基因被用作训练基因。
# 在 training_genes 字段下标注每个 AnnData 中的训练基因，位于 uns 字典中
# 确保数据集中基因顺序的一致性（Tangram 要求每个矩阵中的j列对应相同的基因）。
# 如果一个数据集中某个基因的计数全部为零，则该基因将从训练基因中移除。
# 如果一个基因在两个数据集中都不存在，则该基因将从训练基因中移除。
# 在 pp_adatas 函数中，基因名称被转换为小写以消除不一致的大小写。如果不想这样，可以设置参数 gene_to_lowercase = False
tg.pp_adatas(adata_sc, adata_st, genes=markers)

# Find alignment  寻找比对
# 为 scRNA-seq 谱找到最佳空间对齐，我们使用 map_cells_to_space 函数
# 函数按照 num_epochs 的指定进行迭代映射。我们通常在score达到平台期后中断映射。
# 分数用于衡量映射细胞的基因表达与训练基因的空间数据之间的相似性。
# 默认的映射模式是 mode='cells' ，建议在 GPU 上运行。
# 或者，可以指定 mode='clusters' ，该模式对属于同一簇的单细胞进行平均（通过 cluster_label 传递注释）。这更快，当 scRNAseq 数据和空间数据来自不同标本时，这是我们的选择。
# 如果您希望使用 GPU 运行 Tangram，请设置 device=cuda:0 ，否则使用 device=cpu 
# density_prior 指定每个空间体素内的细胞密度。如果空间体素处于单细胞分辨率（即 MERFISH），请使用 uniform 。默认值 rna_count_based 假定细胞密度与 RNA 分子数量成正比。

ad_map = tg.map_cells_to_space(adata_sc, adata_st,
    mode="cells",
#     mode="clusters",
#     cluster_label='cell_subclass',  # .obs field w cell types
    density_prior='rna_count_based',
    num_epochs=500,
    device='cuda:0',
)

# 映射结果存储在返回的 AnnData 结构中，保存为 ad_map ，结构如下：
# obs 数据框包含单个细胞的元数据。
# var 数据框包含空间数据的元数据。
# uns 字典包含一个数据框，其中包含有关训练基因的各种信息（保存为 train_genes_df ）。

# Cell type maps  细胞类型图谱
# 为了在空间中可视化细胞类型，我们调用 project_cell_annotation 将 annotation 从图谱转移到空间。
# 然后我们可以调用 plot_cell_annotation 进行可视化。
# 您可以通过设置 perc 参数来设定 colormap 的范围，这将有助于去除异常值。

ad_map

tg.project_cell_annotations(ad_map, adata_st, annotation="cell_subclass")
annotation_list = list(pd.unique(adata_sc.obs['cell_subclass']))
tg.plot_cell_annotation_sc(adata_st, annotation_list,perc=0.02)
# 获得映射是否成功的初步判断方法之一是寻找已知的细胞类型模式。为了获得更深入的理解，我们可以使用辅助工具 plot_training_scores ，它提供了四个面板：
tg.plot_training_scores(ad_map, bins=20, alpha=.5)
# 第一个面板是每个训练基因的相似度分数的直方图。
# 第二个面板中，每个点代表一个训练基因，我们可以观察到每个基因的训练分数（y 轴）和在 scRNA-seq 数据中的稀疏性（x 轴）。
# 第三个面板与第二个面板类似，但包含空间数据的基因稀疏性。空间数据通常比单细胞数据更稀疏，这种差异往往是映射质量低下的原因。
# 在最后一个面板中，我们展示了训练分数作为数据集稀疏性差异的函数。对于稀疏性相当的基因，映射的基因表达与空间数据非常相似。然而，如果一个基因在一个数据集（通常是空间数据）中相当稀疏，但在其他数据集中不是，那么映射分数会较低。这发生的原因是 Tangram 由于数据集之间 dropout 数量不一致，无法正确匹配基因模式。

# 尽管上述图表为我们提供了单基因水平的得分概览，但我们仍需了解哪些基因的映射得分较低。
# 这些信息存储在数据框 .uns['train_genes_df'] 中；这是用于构建上述四个图表的数据框。

ad_map.uns['train_genes_df']

# New spatial data via aligned single cells 新的空间数据通过对齐的单个细胞

# 如果映射模式为 'cells' ，我们现在可以使用映射的单个细胞生成"新的空间数据"：这是通过 project_genes 完成的。
# 该函数接受映射（ adata_map ）和相应的单个细胞数据（ adata_sc ）作为输入。结果是一个体素-基因 AnnData ，形式上类似于 adata_st ，但包含来自映射的单个细胞数据的基因表达，而不是 Visium。
# 对于下游分析，我们始终用相应的 ad_ge 替换 adata_st 。
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
ad_ge

# 让我们选择几个映射分数较低的训练基因，尝试理解原因。
genes = ['rragb', 'trim17', 'eno1b']
ad_map.uns['train_genes_df'].loc[genes]
# 为了可视化基因模式，我们使用辅助函数 plot_genes 。该函数接受两个体素-基因 AnnData ：实际空间数据（ adata_measured ）和一个 Tangram 空间预测（ adata_predicted ）。
# 该函数返回两个空间 AnnData 在基因 genes 上的基因表达图。
tg.plot_genes_sc(genes, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)
# 上述图片解释了训练分数较低的原因。一些基因检测到的稀疏度差异很大——通常在空间数据中比在 scRNAseq 中稀疏得多。
# 这是由于像 Visium 这样的技术更容易出现技术性缺失。因此，Tangram 无法为这些基因找到良好的空间对齐，因为基线信号缺失。
# 然而，只要大多数训练基因以高质量测量，我们就可以信任映射并使用 Tangram 预测来校正基因表达。这是一种插补方法，它依赖于与概率模型完全不同的前提。

# 另一种应用是通过检查在空间数据中未被检测到但在单细胞数据中被检测到的基因来发现。
# 它们在使用 pp_adatas 函数进行训练前被移除，但 Tangram 可以预测它们的表达。
genes=['loc102633833', 'gm5700', 'gm8292']
tg.plot_genes_sc(genes, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)
# 迄今为止，我们仅检查了用于对齐数据的基因（训练基因），但映射的单细胞数据 ad_ge 包含了整个转录组。其中包括超过 35k 个测试基因。
(ad_ge.var.is_training == False).sum()


