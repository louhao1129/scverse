# ingest
# 使用ingest映射到reference批次

import scanpy as sc
import pandas as pd
import seaborn as sns
import os
sc.settings.verbosity = 1             # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 3), facecolor='white')
print(os.getcwd())
os.chdir('./scverse/')

# PBMCs 3k(已经处理)
# adata_ref = sc.datasets.pbmc3k_processed() # 下载数据
adata_ref = sc.read_h5ad("./data/pbmc3k_processed.h5ad")
adata = sc.datasets.pbmc68k_reduced()

# To use sc.tl.ingest, the datasets need to be defined on the same variables.
# 求交集
var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names].copy()
adata = adata[:, var_names].copy()

sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
adata_ref

sc.pl.umap(adata_ref, color='louvain')
# 使用 ingest 进行 PBMCs 的映射
sc.tl.ingest(adata, adata_ref, obs='louvain')
# 参数embedding_method：
# 使用adata_ref的umap或pca参数（即adata_ref.obsm['X_pca']和adata_ref.obsm['umap']）去映射得到adata的embedding。

# 参数labeling_method：
# 使用knn去映射得到adata的label（训练数据或参考数据为adata_ref的embedding，测试数据为映射后的adata的embedding，根据相似度赋予adata标签）。

# 参数obs：
# 映射使用的标签键，比如cell_type，louvain等（从adata_ref映射到adata中）。
# 在调用ingest前，需要对adatat_ref运行neighbors()。（neighbors()用于计算观测的neighborhood graph，这是ingest计算前的处理）

adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors'] # fix colors
adata.uns['louvain_colors']

sc.pl.umap(adata, color=['louvain', 'bulk_labels'], wspace=0.5)

df1=adata_ref.obs

adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
df2=adata_concat.obs
## 多了一列信息--batch 区分了数据来源

adata_concat.obs.louvain=adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories )
adata_concat.obs.louvain.cat.categories

adata_concat.uns['louvain_colors'] = adata_ref.uns['louvain_colors']
## 画图展示
sc.pl.umap(adata_concat, color=['batch', 'louvain'])





############
# BBKNN

sc.tl.pca(adata_concat)

sc.external.pp.bbknn(adata_concat, batch_key='batch')

sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['batch', 'louvain'])

############

import scanpy as sc
import pandas as pd
import seaborn as sns
sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')

adata_all = sc.read('F://spe.lesson//k7/pancreas.h5ad' )

adata_all.shape


counts = adata_all.obs.celltype.value_counts()
counts


minority_classes = counts.index[-5:].tolist()        # get the minority classes
adata_all = adata_all[                               # actually subset
    ~adata_all.obs.celltype.isin(minority_classes)]  # ~表示取反
adata_all.obs.celltype.cat.reorder_categories(       # reorder according to abundance
    counts.index[:-5].tolist(), inplace=True)

counts = adata_all.obs.celltype.value_counts()
counts
## 少了五个细胞类型

sc.pp.pca(adata_all)
sc.pp.neighbors(adata_all)
sc.tl.umap(adata_all)
# palette用于选择颜色
sc.pl.umap(adata_all, color=['batch', 'celltype'], palette=sc.pl.palettes.vega_20_scanpy)


sc.external.pp.bbknn(adata_all, batch_key='batch')
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['batch', 'celltype'])
