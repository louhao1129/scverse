# https://scanpy.readthedocs.io/en/stable/tutorials/spatial/integration-scanorama.html
# https://www.nature.com/articles/s41587-019-0113-3

import os
import scanpy as sc
import anndata as an
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanorama

# sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

os.getcwd()
os.chdir('./scverse')
# adata_spatial_anterior = sc.datasets.visium_sge(
#     sample_id="V1_Mouse_Brain_Sagittal_Anterior") # 数据集1下载
# adata_spatial_anterior.write("adata_spatial_anterior.h5ad",compression='gzip')
adata_spatial_anterior=sc.read('adata_spatial_anterior.h5ad')

# adata_spatial_posterior = sc.datasets.visium_sge(
#     sample_id="V1_Mouse_Brain_Sagittal_Posterior") # 数据集2下载
# adata_spatial_posterior.write("adata_spatial_posterior.h5ad",compression='gzip')
adata_spatial_posterior=sc.read('adata_spatial_posterior.h5ad')

adata_spatial_anterior.var_names_make_unique()
adata_spatial_posterior.var_names_make_unique()

adata_spatial_anterior
sc.pp.calculate_qc_metrics(adata_spatial_anterior, inplace=True)
sc.pp.calculate_qc_metrics(adata_spatial_posterior, inplace=True)



for name, adata in [
    ("anterior", adata_spatial_anterior),
    ("posterior", adata_spatial_posterior),]:
    fig, axs = plt.subplots(1, 4, figsize=(12, 3))
    fig.suptitle(f"Covariates for filtering: {name}")

    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.histplot(
        adata.obs["total_counts"][adata.obs["total_counts"] < 20000],
        kde=False,bins=40,ax=axs[1],)
    sns.histplot(adata.obs["n_genes_by_counts"], 
                 kde=False, bins=60, ax=axs[2])
    sns.histplot( adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
        kde=False,bins=60,ax=axs[3],)
        

for adata in [
    adata_spatial_anterior,
    adata_spatial_posterior,
]:
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace=True)


adata_spatial_anterior
# AnnData object with n_obs × n_vars = 2695 × 32285
#     obs: 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes'
#     var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#     uns: 'spatial', 'log1p', 'hvg'
#     obsm: 'spatial'
adata_spatial_posterior
# AnnData object with n_obs × n_vars = 3355 × 32285
#     obs: 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes'
#     var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#     uns: 'spatial', 'log1p', 'hvg'
#     obsm: 'spatial'


adatas = [adata_spatial_anterior, adata_spatial_posterior]
adatas_cor = scanorama.correct_scanpy(adatas, return_dimred=True) # 整合


adatas_cor
# [AnnData object with n_obs × n_vars = 2695 × 32285
#      obs: 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes'
#      var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#      uns: 'spatial', 'log1p', 'hvg'
#      obsm: 'spatial', 'X_scanorama',
#  AnnData object with n_obs × n_vars = 3355 × 32285
#      obs: 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes'
#      var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
#      uns: 'spatial', 'log1p', 'hvg'
#      obsm: 'spatial', 'X_scanorama']
     
    
len(adatas_cor) # 2
adatas_cor[0].obsm['X_scanorama'].shape # (2695, 100)


adata_spatial = sc.concat(
    adatas_cor,label="library_id",
    uns_merge="unique",keys=[
        k
        for d in [
            adatas_cor[0].uns["spatial"],
            adatas_cor[1].uns["spatial"],
        ]
        for k, v in d.items()
    ],
    index_unique="-",
)


adata_spatial



sc.pp.neighbors(adata_spatial, use_rep="X_scanorama")
sc.tl.umap(adata_spatial)
sc.tl.leiden(adata_spatial, key_added="clusters")

sc.pl.umap(
    adata_spatial, color=["clusters", "library_id"], palette=sc.pl.palettes.default_20)


clusters_colors = dict(
    zip([str(i) for i in range(18)], adata_spatial.uns["clusters_colors"]))
    
    


fig, axs = plt.subplots(1, 2, figsize=(15, 10))

for i, library in enumerate(
    ["V1_Mouse_Brain_Sagittal_Anterior", "V1_Mouse_Brain_Sagittal_Posterior"]
):
    ad = adata_spatial[adata_spatial.obs.library_id == library, :].copy()
    sc.pl.spatial(
        ad,img_key="hires", library_id=library,
        color="clusters",size=1.5,palette=[
            v
            for k, v in clusters_colors.items()
            if k in ad.obs.clusters.unique().tolist()],
        legend_loc=None,
        show=False,
        ax=axs[i], )

plt.tight_layout()


os.getcwd()
counts = pd.read_csv("geo\\GSE115746_cells_exon_counts.csv", index_col=0).T
meta = pd.read_csv("geo\\GSE115746_complete_metadata_28706-cells.csv", index_col="sample_name")
meta = meta.loc[counts.index]

annot = sc.queries.biomart_annotations(
        "mmusculus",
        ["mgi_symbol", "ensembl_gene_id"],
    ).set_index("mgi_symbol")
    
annot = annot[annot.index.isin(counts.columns)]
counts = counts.rename(columns=dict(zip(annot.index, annot["ensembl_gene_id"])))
adata_cortex = an.AnnData(counts, obs=meta)

sc.pp.normalize_total(adata_cortex, inplace=True)
sc.pp.log1p(adata_cortex)
adata_cortex.write_h5ad("geo\\adata_processed.h5ad")

adata_spatial_anterior.var.set_index("gene_ids", inplace=True) # 将 "gene_ids" 列设置为索引（index）
adata_spatial_posterior.var.set_index("gene_ids", inplace=True)

adata_anterior_subset = adata_spatial_anterior[
    adata_spatial_anterior.obsm["spatial"][:, 1] < 6000, :]
adata_posterior_subset = adata_spatial_posterior[
 (adata_spatial_posterior.obsm["spatial"][:, 1] < 4000)
 & (adata_spatial_posterior.obsm["spatial"][:, 0] < 6000),:,]

adatas_anterior = [adata_cortex, adata_anterior_subset]
adatas_posterior = [adata_cortex, adata_posterior_subset]

adatas_cor_anterior = scanorama.correct_scanpy(adatas_anterior, return_dimred=True)
adatas_cor_posterior = scanorama.correct_scanpy(adatas_posterior, return_dimred=True)


adata_cortex_anterior = sc.concat(
    adatas_cor_anterior,label="dataset",
    keys=["smart-seq", "visium"],join="outer",
    uns_merge="first",)
adata_cortex_posterior = sc.concat(
    adatas_cor_posterior,label="dataset",
    keys=["smart-seq", "visium"],join="outer",
    uns_merge="first",)


# see paper：https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8

from sklearn.metrics.pairwise import cosine_distances

distances_anterior = 1 - cosine_distances(
    adata_cortex_anterior[adata_cortex_anterior.obs.dataset == "smart-seq"].obsm[
        "X_scanorama"
    ],
    adata_cortex_anterior[adata_cortex_anterior.obs.dataset == "visium"].obsm[
        "X_scanorama"
    ],
)
distances_posterior = 1 - cosine_distances(
    adata_cortex_posterior[adata_cortex_posterior.obs.dataset == "smart-seq"].obsm[
        "X_scanorama"
    ],
    adata_cortex_posterior[adata_cortex_posterior.obs.dataset == "visium"].obsm[
        "X_scanorama"
    ],
)


def label_transfer(dist, labels):
    lab = pd.get_dummies(labels).to_numpy().T
    class_prob = lab @ dist
    norm = np.linalg.norm(class_prob, 2, axis=0)
    class_prob = class_prob / norm
    class_prob = (class_prob.T - class_prob.min(1)) / class_prob.ptp(1)
    return class_prob

class_prob_anterior = label_transfer(distances_anterior, adata_cortex.obs.cell_subclass)
class_prob_posterior = label_transfer(
    distances_posterior, adata_cortex.obs.cell_subclass)

cp_anterior_df = pd.DataFrame(
    class_prob_anterior,
    columns=sorted(adata_cortex.obs["cell_subclass"].cat.categories),)
cp_posterior_df = pd.DataFrame(
    class_prob_posterior,
    columns=sorted(adata_cortex.obs["cell_subclass"].cat.categories),)

cp_anterior_df.index = adata_anterior_subset.obs.index
cp_posterior_df.index = adata_posterior_subset.obs.index


adata_anterior_subset_transfer = adata_anterior_subset.copy()
adata_anterior_subset_transfer.obs = pd.concat(
    [adata_anterior_subset.obs, cp_anterior_df], axis=1
)

adata_posterior_subset_transfer = adata_posterior_subset.copy()
adata_posterior_subset_transfer.obs = pd.concat(
    [adata_posterior_subset.obs, cp_posterior_df], axis=1
)

sc.pl.spatial(
    adata_anterior_subset_transfer,
    img_key="hires",
    color=["L2/3 IT", "L4", "L5 PT", "L6 CT"],
    size=1.5,
)
sc.pl.spatial(
    adata_posterior_subset_transfer,
    img_key="hires",
    color=["L2/3 IT", "L4", "L5 PT", "L6 CT"],
    size=1.5,
)

sc.pl.spatial(adata_anterior_subset_transfer,
     img_key="hires", color=["Oligo", "Astro"], size=1.5)


sc.pl.spatial(adata_posterior_subset_transfer, 
    img_key="hires", color=["Oligo", "Astro"], size=1.5)

adata_posterior_subset_transfer.write_h5ad("F:\\spe.lesson\\k8\\geo\\adata_posterior_subset_transfer.h5ad")
adata_anterior_subset_transfer.write_h5ad("F:\\spe.lesson\\k8\\geo\\adata_anterior_subset_transfer.h5ad")

