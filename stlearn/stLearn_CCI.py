import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir('./scverse/' )

# 数据加载与预处理
data_dir = "./data/stlearn/BCBA"

data=st.Read10X(data_dir)
data.var_names_make_unique()
st.add.image(
    adata=data,imgpath=data_dir+"/spatial/tissue_hires_image.png",
    library_id="V1_Breast_Cancer_Block_A_Section_1",visium=True
)


st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data) # NOTE: no log1p


# Adding the label transfer results
spot_mixtures = pd.read_csv(data_dir+'/label_transfer_bc.csv', index_col=0, sep='\t')
labels = spot_mixtures.loc[:,'predicted.id'].values.astype(str)
spot_mixtures = spot_mixtures.drop(['predicted.id','prediction.score.max'], axis=1)
spot_mixtures.columns = [col.replace('prediction.score.', '') for col in spot_mixtures.columns]
# Note the format! #
print(labels)

print(spot_mixtures)


# Check is in correct order
print('Spot mixture order correct?"',np.all(spot_mixtures.index.values==data.obs_names.values)) # Check is in correct order


# NOTE: using the same key in data.obs & data.uns
data.obs['cell_type'] = labels # Adding the dominant cell type labels per spot
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures # Adding the cell type scores


st.pl.cluster_plot(data, use_label='cell_type')

# 运行配体-受体分析
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
lrs

# Running the analysis #
st.tl.cci.run(data, lrs,
					min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots                  
					distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode                  
					n_pairs=10000, # Number of random pairs to generate; recommend ~10,000                  
					n_cpus=32) # Number of CPUs for parallel. If None, detects & use all available.                  



# A dataframe detailing the LR pairs ranked by number of significant spots.
lr_info = data.uns['lr_summary']
print('\n', lr_info)

# P 值矫正
st.tl.cci.adj_pvals(data, correct_axis='spot',pval_adj_cutoff=0.05, adj_method='fdr_bh')

# 可视化 LRs 在显著位点上的整体排名
st.pl.lr_summary(data, n_top=500)
st.pl.lr_summary(data, n_top=50, figsize=(10,3))

# 诊断
# A key aspect of the LR analysis is to control for LR expression level and frequency when calling significant hotspots.
#  our diagnostic plots should show next to no correlation between the hotspots of the LR pairs and the expression level and frequency of expression
# if not, could indicate a larger number of permutations is required.
st.pl.lr_diagnostics(data, figsize=(10,2.5))

# Left plot: Relationship between LR expression level (non-zero spots average median expression of genes in the LR pair) and the ranking of the LR.
# Right plot: Relationship between LR expression frequency (average proportion of zero spots for each gene in the LR pair) and the ranking of the LR.

#In this case, there is a weak correlation between the LR expression frequency and number of significant spots, 
# indicating the n_pairs parameter should be set higher to create more accurate background distributions (10,000 pairs was used in the case of the paper version of the above).
st.pl.lr_n_spots(data, n_top=50, figsize=(11, 3), max_text=100)
st.pl.lr_n_spots(data, n_top=500, figsize=(11, 3), max_text=100)
# The above boxplots show the number of spots with ligand-receptor scores for each LR on the y-axis, with the LR ranking on the x-axis. 
# The bar area is stratified by spots which were significant (green) and non-significant (blue).
# While the trend isn’t as quantitative with this plot compared with the scatter plot, 
# there does appear to be some correlation with more highly frequent LR pairs and LR ranking; again indicating the n_pairs parameter above should be set higher.



# Biologic Plots (Optional)

## Running the GO enrichment analysis ##
# r_path = "/home/louhao/.conda/envs/rl/bin/R"
# st.tl.cci.run_lr_go(data, r_path)

# st.pl.lr_go(data, lr_text_fp={'weight': 'bold', 'size': 10}, rot=15,
#                figsize=(12,3.65), n_top=15, show=False)

# 总体而言，我们观察到一些较强的生物学富集，表明存在一些由排名靠前的 LR pairs 介导的潜在通路。

# LR Statistics Visualisations
best_lr = data.uns['lr_summary'].index.values[0] # Just choosing one of the top from lr_summary
stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
fig, axes = plt.subplots(ncols=len(stats), figsize=(16,6))
for i, stat in enumerate(stats):
    st.pl.lr_result_plot(data, use_result=stat, use_lr=best_lr, show_color_bar=False, ax=axes[i])
    axes[i].set_title(f'{best_lr} {stat}')

# This shows the scores -log10(p_adjs) for all spots, then the scores subseted to significant spots
fig, axes = plt.subplots(ncols=2, figsize=(8,6))
st.pl.lr_result_plot(data, use_result='-log10(p_adjs)', use_lr=best_lr, show_color_bar=False, ax=axes[0])
st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr, show_color_bar=False, ax=axes[1])
axes[0].set_title(f'{best_lr} -log10(p_adjs)')
axes[1].set_title(f'{best_lr} lr_sig_scores')
# LR Interpretation Visualisations
# These visualisations are meant to help interpret the directionality of the cross-talk.

# Binary LR coexpression plot for all spots
st.pl.lr_plot(data, best_lr, inner_size_prop=0.1, outer_mode='binary', pt_scale=5,
              use_label=None, show_image=True,
              sig_spots=False)
# Binary LR coexpression plot for significant spots
st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode='binary', pt_scale=20,
              use_label=None, show_image=True,
              sig_spots=True)
# Continuous LR coexpression for significant spots
# The receptor is in green, the ligand is in red. The inner-point is the receptor, the outter point is the ligand.
# Help to see where and how heavily expressed ligands/receptors are.
# Idea is receptor is on the cell surface, & ligand permeates out from the cell surface.
# All spots #
st.pl.lr_plot(data, best_lr,
              inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,
              outer_mode='continuous', pt_scale=60,
              use_label=None, show_image=True,
              sig_spots=False)
# Only significant spots #
st.pl.lr_plot(data, best_lr,
              inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,
              outer_mode='continuous', pt_scale=60,
              use_label=None, show_image=True,
              sig_spots=True)
# .......

# 预测显著细胞间相互作用

# 诊断图：检查交互作用和细胞类型频率相关性

# 细胞间相互作用可视化