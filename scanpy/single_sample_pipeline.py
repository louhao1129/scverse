# Core scverse libraries
import scanpy as sc
import anndata as ad
import os
sc.settings.set_figure_params(dpi=150, facecolor="white")
os.getcwd()
os.chdir("./scverse/")

# 数据下载地址：https://figshare.com/articles/dataset/NeurIPS_2021_Benchmark_dataset/22716739
#设置文件的名字
samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}
adatas = {}

for sample_id, filename in samples.items():
    path = f"./data/{filename}"
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique() # To make them unique
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
print(adata.obs["sample"].value_counts())
adata

# 质量控制
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

## QC 指标的小提琴图
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# 根据质控指标图，现在可以通过设置手动或自动阈值来移除表达过多线粒体基因或总计数过多的细胞。
# 然而，有时看似较差的质控指标可能是由真实生物学现象驱动的，
# 因此我们建议从非常宽松的过滤策略开始，并在稍后阶段重新审视。
# 因此，我们现在仅过滤表达少于 100 个基因的细胞以及在少于 3 个细胞中检测到的基因。


# 此外，需要注意的是，对于包含多个批次的数据集，应针对每个样本单独进行质量控制，
# 因为批次之间的质量控制阈值可能会有很大差异。
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# 双细胞检测

sc.pp.scrublet(adata, batch_key="sample")
# 在obs中生成新的列，表示每个细胞的双细胞分数和预测结果
# 我们可以通过过滤被标记为双细胞的细胞，或在完成聚类后过滤具有高双细胞分数的簇来去除双细胞。
# n0 = adata.shape[0]
# adata = adata[adata.obs['predicted_doublet']==False, :].copy()
# n1 = adata.shape[0]
# print(f'Cells retained after scrublet: {n1}, {n0-n1} removed.')

# 归一化
# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata) # Seurat use target_sum=1e4
# Logarithmize the data
sc.pp.log1p(adata)

# Feature selection
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata)
# 高变基因检测的结果被存储为注释信息在 .var["highly_variable"] 中
# 并由 PCA 自动检测，因此 sc.pp.neighbors 以及后续的流形/图工具也会自动检测。
# 在这种情况下，实际上执行以下过滤步骤也是不必要的。

# 降维
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# UMAP可视化
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)
# 尽管本教程考虑的数据包含两个不同的样本，
# 我们仅观察到轻微的批次效应，并可以继续进行数据的聚类和注释。

# 如果你的 UMAP 图展示出批次效应，整合样本并执行批次校正/整合可能是有益的
# 批次整合方法：scanorama、harmony、bbknn、scvi-tools等

# 聚类
# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])

# 重新评估质量控制与细胞过滤
# 我们将通过使用 UMAP 可视化不同的 QC 指标来重新评估我们的过滤策略
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)
sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

# 手动细胞类型注释 
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)
# Marker 基因
# 我们为这个数据集中预期的主要细胞类型定义一组标记基因
marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    # Note HBM and GYPA are negative markers
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    # Note IGHD and IGHM are negative markers
    "B cells": [
        "MS4A1",
        "ITGB1",
        "COL4A4",
        "PRDM1",
        "IRF4",
        "PAX5",
        "BCL11A",
        "BLK",
        "IGHD",
        "IGHM",
    ],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    # Note PAX5 is a negative marker
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")

adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Lymphocytes",
        "1": "Monocytes",
        "2": "Erythroid",
        "3": "B Cells",
    }
)

# 差异表达基因作为marker基因
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
# 在散点图上可视化表达差异最大的前 5 个基因
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5 # min-max scale
)
# https://github.com/scverse/scanpy/issues/1757 ： sc.pl.dotplot() standard_scale='var': should scaling be changed? 
# 我们可以利用这些基因来确定我们正在观察的细胞类型。例如，簇 7 表达 NKG7 和 GNLY，表明这些是 NK 细胞
# 差异表达基因可以以方便的格式提取
sc.get.rank_genes_groups_df(adata, group="7").head(5)

dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
)
# 您可能已经注意到，这里得到的 p 值非常低。这是由于所进行的统计测试将每个细胞视为一个独立样本。
# 为了采取更保守的方法，您可以考虑按样本（例如 sc.get.aggregate(adata, by=["sample", "cell_type"], func="sum", layer="counts") ）
# 进行“伪合并”您的数据，并使用更强大的差异表达工具，如 pydeseq2 。

adata.write("./data/pbmc3k_preprocess.h5ad")
