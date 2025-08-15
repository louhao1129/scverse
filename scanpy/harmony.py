# https://github.com/slowkow/harmonypy

import scanpy as sc
import pandas as pd
import seaborn as sns
import os
sc.settings.verbosity = 1             # verbosity errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(3, 3), facecolor='white')

print(os.getcwd())
os.chdir('./k8')

data_s10=sc.read_10x_mtx('multiple/BC2')
data_s21=sc.read_10x_mtx('multiple/BC10')
data_s2=sc.read_10x_mtx('multiple/BC21')

data_s10.obs["orig"]='s10'
data_s21.obs["orig"]='s21'
data_s2.obs["orig"]='s2'

adatas = sc.AnnData.concatenate(data_s10,data_s21,data_s2)
adatas.obs

sc.pp.filter_cells(adatas, min_genes=200)
sc.pp.filter_genes(adatas, min_cells=3)

adatas.var['mt'] = adatas.var_names.str.startswith('MT-')
# adatas.var['mt']
sc.pp.calculate_qc_metrics(adatas, qc_vars=['mt'], 
    percent_top=None, log1p=False, inplace=True)
adatas.obs

sc.pl.violin(adatas, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

adatas = adatas[adatas.obs.n_genes_by_counts > 500, :]
adatas = adatas[adatas.obs.pct_counts_mt < 20, :]

sc.pp.normalize_total(adatas, target_sum=1e4)
sc.pp.log1p(adatas)
sc.pp.highly_variable_genes(adatas, min_mean=0.0125, max_mean=3, min_disp=0.5)

adatas.raw = adatas
adatas = adatas[:, adatas.var.highly_variable]

sc.pp.regress_out(adatas, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(adatas, max_value=10)
sc.pp.pca(adatas)
sc.pp.neighbors(adatas)
sc.tl.umap(adatas)
sc.pl.umap(adatas, color=['batch', 'orig'], palette=sc.pl.palettes.vega_20_scanpy)

import scanpy.external as sce
sce.pp.harmony_integrate(adatas, 'orig') # harmony
sc.pp.neighbors(adatas, use_rep="X_pca_harmony")
sc.tl.umap(adatas,init_pos='X_pca_harmony')
sc.pl.umap(adatas, color=['batch', 'orig'], legend_fontsize=8)

adatas.obsm

sc.tl.leiden(adatas )
sc.pl.umap(adatas, color=['leiden', 'CD3D' ])

# bbknn
sc.external.pp.bbknn(adatas, batch_key='batch')
sc.tl.umap(adatas)
sc.pl.umap(adatas, color=['batch', 'orig'])
adatas.obsm
sc.tl.leiden(adatas )

sc.pl.umap(adatas, color=['leiden', 'CD3D', 'NKG7'])


sc.tl.rank_genes_groups(adatas, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adatas, n_genes=25, sharey=False)
pd.DataFrame(adatas.uns['rank_genes_groups']['names']).head(5)


new_cluster_names = [
    '0', '1', '2: CD14+Monocytes','3 : epi', '4:endo',
    '5: NK', '6;T', '7;T', '8: Megakaryocytes','9',
    '10','11','12','13','14','15','16','17','18',
    '19','20','21','22','23','24']
adatas.rename_categories('leiden', new_cluster_names)


sc.pl.umap(adatas, color='leiden', legend_loc='on data', title='', frameon=False)

sc.pl.rank_genes_groups_tracksplot(adatas, n_genes=3)

