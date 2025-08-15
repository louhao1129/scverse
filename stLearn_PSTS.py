# Pseudotime Spatial Trajectory Inference
# 在本教程中，我们使用空间信息以及基因表达谱来执行空间轨迹推断，以探索导管内癌（DCIS）向浸润性导管癌（IDC）的进展过程
import stlearn as st
import scanpy as sc
import pathlib as pathlib
import os
import warnings
warnings.filterwarnings("ignore") # Ignore all warnings
st.settings.set_figure_params(dpi=120)


print(os.getcwd())
os.chdir("./scverse/")
# 1. 准备工作

# 读取和预处理数据

# Read raw data
st.settings.datasetdir = pathlib.Path.cwd() / "data" / "stlearn"
# data = st.datasets.visium_sge(sample_id="V1_Breast_Cancer_Block_A_Section_1")
data = sc.read_h5ad("./data/stlearn/V1_Breast_Cancer_Block_A_Section_1.h5ad")

# stLearn 也像 scanpy 一样使用 AnnData 作为核心对象。但是，我们以不同的方式存储空间数据。
# 如果你通过 scanpy 读取数据，需要将其转换为 stLearn AnnData 格式，即使用data = st.convert_scanpy(data)
# 但是，如果你通过 stLearn 读取数据，你可以使用 scanpy 的几乎所有函数。
data = st.convert_scanpy(data)

# Save raw_count
data.layers["raw_count"] = data.X
# Preprocessing
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
# Keep raw data
data.raw = data
st.pp.scale(data)

# 聚类数据
# Run PCA
st.em.run_pca(data, n_comps=50, random_state=0)
# Tiling image, 图像分割
st.pp.tiling(data, out_path="tiling", crop_size=40)
# Using Deep Learning to extract feature， 深度学习识别图像特征
st.pp.extract_feature(data)
# Apply stSME spatial-PCA option， SME标准化数据
# stSME 是一种新型归一化方法，它专为空间转录组学数据设计，并利用组织的空间位置、形态以及基因表达。
# 一种基于空间图的神经网络插补方法（stSME），用于校正技术噪声/缺失值并增加 ST 数据覆盖范围
st.spatial.morphology.adjust(data, use_data="X_pca", radius=50, method="mean")
st.pp.neighbors(data, n_neighbors=25, use_rep='X_pca_morphology', random_state=0)
st.tl.clustering.louvain(data, random_state=0)

st.pl.cluster_plot(data, use_label="louvain", image_alpha=1, size=7)

st.add.annotation(data, label_list=['Fatty tissue,immune/lymphoid 1 MALAT1+',
                                    'Invasive cancer,fibrous tissue 1 CXCL14+',
                                    'Invasive cancer,fibrous tissue 2 CRISP3+',
                                    'Invasive cancer,fibrous tissue, fatty tissue',
                                    'Fatty tissue,immune/lymphoid 2 IGKC+',
                                    'Fibrous tissue',
                                    'Invasive cancer,fibrous tissue (DCIS)',
                                    'Fatty tissue, Fibrous tissue',
                                    'Invasive cancer,immune/lymphoid (IDC)',
                                    'Invasive cancer,fatty tissue 3 MUC5B+',
                                    'Fatty tissue'],
                  use_label="louvain")
st.pl.cluster_plot(data, use_label="louvain_anno", image_alpha=1, size=7)

# 2. 空间轨迹推断

# Choosing root  选择根节点
# 3733 是我们选择的根节点的索引。它位于 DCIS 聚类（6）中。
# 我们建议根节点应该位于 UMAP 空间中聚类的末端/起始端。你可以在 UMAP 中找到聚类的最小/最大点作为根节点

data.uns["iroot"] = st.spatial.trajectory.set_root(data, use_label="louvain", cluster="6", use_raw=True)
st.spatial.trajectory.pseudotime(data, eps=50, use_rep="X_pca", use_label="louvain")

# 空间轨迹推断 - 全局层面
# 我们运行全局层面的伪时间-空间（PSTS）方法来重建聚类 6（DCIS）和 8（IDC 病变）之间的空间轨迹。
st.spatial.trajectory.pseudotimespace_global(data, use_label="louvain", list_clusters=["6", "8"])

st.pl.cluster_plot(data, use_label="louvain", show_trajectories=True, list_clusters=["6", "8"], show_subcluster=True)
st.pl.trajectory.tree_plot(data)

# 检测轨迹相关marker基因
# 根据空间轨迹/树状图，我们可以看到 2 个分支从子簇 6 和 15 开始。然后我们运行函数来检测与 PSTS 值高度相关的基因。
st.spatial.trajectory.detect_transition_markers_clades(data, clade=6, use_raw_count=False, cutoff_spearman=0.4)
# Detecting the transition markers of clade_6...
# Transition markers result is stored in adata.uns['clade_6']
st.spatial.trajectory.detect_transition_markers_clades(data, clade=15, use_raw_count=False, cutoff_spearman=0.4)
st.spatial.trajectory.detect_transition_markers_clades(data, clade=21, use_raw_count=False, cutoff_spearman=0.4)

# 剔除核糖体基因
data.uns['clade_6'] = data.uns['clade_6'][data.uns['clade_6']['gene'].map(lambda x: "RPL" not in x)]
data.uns['clade_15'] = data.uns['clade_15'][data.uns['clade_15']['gene'].map(lambda x: "RPL" not in x)]
data.uns['clade_21'] = data.uns['clade_21'][data.uns['clade_21']['gene'].map(lambda x: "RPL" not in x)]
data.uns['clade_6'] = data.uns['clade_6'][data.uns['clade_6']['gene'].map(lambda x: "RPS" not in x)]
data.uns['clade_15'] = data.uns['clade_15'][data.uns['clade_15']['gene'].map(lambda x: "RPS" not in x)]
data.uns['clade_21'] = data.uns['clade_21'][data.uns['clade_21']['gene'].map(lambda x: "RPS" not in x)]

# 在轨迹相关marker基因图上，左侧（红色）的基因与空间轨迹呈负相关，而右侧（蓝色）的基因与空间轨迹呈正相关。
# top marker 基因可视化
st.pl.trajectory.transition_markers_plot(data, top_genes=30, trajectory="clade_6")
st.pl.trajectory.transition_markers_plot(data, top_genes=30, trajectory="clade_15")
st.pl.trajectory.transition_markers_plot(data, top_genes=30, trajectory="clade_21")

# 还提供了一个函数，用于比较两个类群之间的转换标记。
st.spatial.trajectory.compare_transitions(data, trajectories=["clade_15", "clade_21"])
st.pl.trajectory.DE_transition_plot(data)