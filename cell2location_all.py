import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cell2location.utils import select_slide
from cell2location.plt import plot_spatial
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

print(os.getcwd())
# os.chdir("./scverse/")
# 注意不要把这个脚本命令为cell2location.py，否则可能出现cell2location模块的导入错误
# 定义结果存储位置
results_folder = './results/lymph_nodes_analysis/'
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# 加载 Visium 和 scRNA-seq 参考数据

# adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata_vis = sc.read_visium("./data/V1_Human_Lymph_Node")
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

# 这里我们将基因重命名为 ENSEMBL ID，以便在单细胞数据和空间数据之间进行正确匹配
# 因此你可以忽略 scanpy 建议调用 .var_names_make_unique 。
adata_vis.var['SYMBOL'] = adata_vis.var_names 
adata_vis.var.set_index('gene_ids', drop=True, inplace=True) # set_index修改行名

# 仍然可以使用标准的 scanpy 函数按名称绘制基因表达，如下所示
sc.pl.spatial(adata_vis, color='PTPRC', gene_symbols='SYMBOL')

# 线粒体编码基因（基因名称以前缀 mt-或 MT-开头）与空间映射无关，因为它们的表达在单细胞和核数据中代表技术伪影，而不是线粒体的生物丰度。
# 然而，这些基因占每个位置 mRNA 的 15-40%。因此，为了避免映射伪影，我们强烈建议去除线粒体基因。

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# Read data
# adata_ref = sc.read(
#     f'./data/sc.h5ad',
#     backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad'
# )
adata_ref = sc.read_h5ad("./data/sc.h5ad")
adata_ref.var['SYMBOL'] = adata_ref.var.index
# rename 'GeneID-2' as necessary for your data
adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)
# delete unnecessary raw slot (to be removed in a future version of the tutorial)
del adata_ref.raw

# 在进行参考细胞类型特征估计之前，我们建议执行非常宽松的基因选择。我们更倾向于标准的高变基因选择，
# 因为我们的方法保留了稀有基因的标记，同时去除了大多数信息量不大的基因。
# 默认参数 cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12 是一个良好的起点，
# 但是，可以增加阈值以排除更多基因。为了保留稀有细胞类型的标记基因，我们建议低 cell_count_cutoff=5 ，
# 但是， cell_percentage_cutoff2 和 nonz_mean_cutoff 可以增加,以选择 8k-16k 个基因。
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# 在这个二维直方图中，橙色矩形突出了基于基因表达细胞数量（Y 轴）和检测到该基因的细胞平均 RNA 计数（X 轴）组合排除的基因。
# 在这种情况下，下载的数据集已经使用这种方法进行了过滤，因此橙色矩形下没有密度（未来教程版本将进行更改）。
# filter the object
adata_ref = adata_ref[:, selected].copy()



# 估计参考细胞类型特征（NB 回归）
# 这些特征是通过使用负二项回归模型从 scRNA-seq 数据中估计的，考虑了批次效应
# 首先，为回归模型准备 anndata 对象
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='Subset',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['Method']
                       )

# create the regression model

mod = RegressionModel(adata_ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()

# 训练模型来估计参考细胞类型特征
# 请注意，为了在您的数据上收敛（=获得损失稳定），您可能需要增加 max_epochs=250
# 还请注意，在这里我们使用 batch_size=2500 ，它比 scvi-tools 默认值大得多，并在数据中的所有细胞上（ train_size=1 ）进行训练——这两个参数都是默认值。
mod.train(max_epochs=250)

# 确定模型是否需要更多训练
# 此处我们绘制了训练过程中的 ELBO 损失历史，从图中去除了前 20 个 epoch。
# 该图应呈现下降趋势，并在训练结束时趋于平稳。如果它仍在下降，请增加 max_epochs
mod.plot_history(20)



# 将训练好的模型用于后验推断，以估计参考细胞类型签名的后验分布
# 从训练好的模型中导出后验样本，这些样本用于表示参数的不确定性
# num_samples：指定要抽样的后验样本的数量
# batch_size：控制每次推断所用的数据批次的大小
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

# 可以直接计算后验分布的 5%、50%和 95%分位数，而不是从分布中抽取 1000 个样本（或任何其他分位数）。
# 这加快了在大型数据集上的应用，并减少了内存需求——然而，后验均值和标准差不能通过这种方式计算。
# adata_ref = mod.export_posterior(
#     adata_ref, use_quantiles=True,
#     # choose quantiles
#     add_to_obsm=["q05","q50", "q95", "q0001"],
#     sample_kwargs={'batch_size': 2500}
# )

# 检查质控图
# 1.评估重建精度以判断是否存在推理问题。这个二维直方图应该有大部分观测值沿着一个有噪声的对角线。
# 2.由于批次效应，估计的表达特征与每个簇中的平均表达不同。对于不受批次效应影响的 scRNA-seq 数据集（本数据集就是这样），可以使用簇平均表达代替使用模型估计特征。当这个图与对角线图非常不同时（例如 Y 轴值非常低，到处都有密度），这表明特征估计存在问题。
mod.plot_QC()

# 模型和输出 h5ad 文件可以按如下方式加载
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)


# 提取参考细胞类型特征作为 pd.DataFrame
# 所有负二项回归模型的参数都被导出到参考 anndata 对象中，但对于空间映射，我们只需要每个细胞类型中每个基因的估计表达量。
# 在这里我们从标准输出中提取这些信息
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# Cell2location：空间映射
# 寻找共享基因并准备 anndata。对 anndata 和参考signature 取子集

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# 要使用 cell2location 空间映射模型，您需要指定 2 个超参数（ N_cells_per_location 和 detection_alpha ）
# ——有关设置这些超参数及其影响的详细说明，请参阅流程图和注释。

# 选择超参数 N_cells_per_location
# 将预期的细胞丰度 N_cells_per_location 适应到每个组织中是有用的。此值可以从配对的病理学图像中估计
# 如上注释所述。将本教程中显示的值（ N_cells_per_location=30 ）更改为您组织中观察到的值。

# 选择超参数 detection_alpha

# 为了提高在 RNA 检测敏感性在载玻片/批次中具有较大技术变异的数据集上的准确性和灵敏度——你需要放宽每个位置的标准化正则化（使用 detection_alpha=20 ）。
# 当您观察到每个位置的 RNA 总计数的空间分布与基于组织学检查的预期细胞数量不匹配时，表明您的样本中存在 RNA 检测敏感性的高技术变异。
# 我们最初选择高正则化（ detection_alpha=200 ）作为默认值，因为我们在论文中使用的鼠脑和人类淋巴结数据集技术效应较低，使用高正则化强度可以提高每个位置的总估计细胞丰度与从组织学定量计数的细胞核数量之间的稳定性（cell2location 论文中的图 S8F）。
# 然而，在许多合作中，我们发现人类组织上的 Visium 实验受到技术效应的影响。这促使我们设定新的默认值 detection_alpha=20 ，
# 并建议在您的数据上测试这两种设置（ detection_alpha=20 和 detection_alpha=200 ）。
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()
# 训练 cell2location
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training'])

# 导出细胞丰度后验分布估计值并保存结果
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

# 模型和输出 h5ad 文件可以按如下方式加载
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# 评估映射质量。检查重建精度以确定是否存在映射问题。该图应大致呈对角线，显著的偏差将表明需要调查的问题。
mod.plot_QC()
# 在进行多个空间批次整合以及处理在切片中具有显著差异的 RNA 检测数据集（这些差异无法用组织学中的高细胞密度来解释）时，评估 cell2location 是否已对这些效应进行标准化非常重要。
# 您期望在不同批次中观察到相似的细胞总丰度，但具有不同的 RNA 检测灵敏度（均由 cell2location 估计）。
# 您期望细胞总丰度能反映组织学中的高细胞密度。
fig = mod.plot_spatial_QC_across_batches()


# 在空间坐标中可视化细胞丰度
# 我们使用后验分布的 5%分位数，代表模型对细胞丰度值有高度置信度（即“至少存在此数量”）。
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# 选择一个切片
slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')
# 在空间坐标系中画图
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(slide, cmap='magma',
                  # show first 8 cell types
                  color=['B_Cycling', 'B_GC_LZ', 'T_CD4+_TfH_GC', 'FDC',
                         'B_naive', 'T_CD4+_naive', 'B_plasma', 'Endo'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )

# 在一张图张展示多个细胞类型
# select up to 6 clusters
clust_labels = ['T_CD4+_naive', 'B_naive', 'FDC']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )

# 下游分析 
# 通过 Leiden 聚类识别离散组织区域
# 我们通过使用 cell2location 估计的细胞丰度对位置进行聚类，来识别细胞组成不同的组织区域。
# 我们通过使用每种细胞类型的估计细胞丰度对 Visium 位点进行聚类来识别组织区域。
# 我们构建了一个表示估计细胞丰度中位置相似性的 K 近邻（KNN）图，然后应用 Leiden 聚类。
# KNN 近邻的数量应根据数据集的大小和解剖学定义区域的大小进行调整（例如，海马区域相对于大脑的大小较小，因此可能被较大的 n_neighbors 掩盖）。
# 这可以通过对一系列 KNN 近邻和 Leiden 聚类分辨率进行操作，直到获得与组织解剖结构匹配的聚类为止。
# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)

# 聚类是在所有 Visium 切片/批次中联合完成的，因此区域身份可以直接进行比较。
# 当多个批次之间存在强技术效应（这个例子并非如此）时， sc.external.pp.bbknn 原则上可用于在 KNN 构建过程中考虑这些效应。

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=1.1)
# 得到的聚类结果保存在 adata_vis.obs['region_cluster'] 中
# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")

# 我们可以使用位置组成相似性图来构建所有切片/Visium 批次的联合集成 UMAP 表示。
# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_vis, min_dist = 0.3, spread = 1)

# show regions in UMAP coordinates
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_vis, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20)
    sc.pl.umap(adata_vis, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20)

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(adata_vis, color=['region_cluster'],
                  size=1.3, img_key='hires', alpha=0.5)

# 使用矩阵分解（NMF）识别细胞区室/组织区域

# 估计空间数据中每个基因的细胞类型特异性表达（NCEM 所需）

# 高级用法

