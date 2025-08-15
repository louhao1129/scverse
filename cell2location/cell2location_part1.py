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
