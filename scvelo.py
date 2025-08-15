import scanpy as sc
import scvelo as scv
import pandas as pd
import os

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

os.getcwd()
os.chdir("./scverse/")
# adata = scv.datasets.pancreas() # 下载数据
adata = sc.read_h5ad("./data/velo/pancreas.h5ad")
adata
scv.pl.proportions(adata)

# 预处理数据
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
sc.pp.log1p(adata)

# 此外，我们需要在 PCA 空间中计算最近邻的一阶和二阶矩（均值和非中心化方差），
# 这些内容总结在 scv.pp.moments 中，该函数内部计算 scv.pp.pca 和 scv.pp.neighbors 。一阶矩对于确定性速度估计是必需的，而随机估计也需要二阶矩。
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# 可以进行进一步预处理（如批次效应校正）以去除不需要的变异来源。有关最佳实践，请参阅详细信息。
# 注意，任何额外的预处理步骤仅影响 X，并且不应用于剪接/未剪接计数。

# 估计 RNA 速度
# 速度是基因表达空间中的矢量，表示单个细胞移动的方向和速度。速度通过建模剪接动力学中的转录动态获得，
# 可以是随机（默认）或确定性（通过设置 mode='deterministic' ）。
# 对于每个基因，拟合了前体成熟（未剪接）和成熟（剪接）mRNA 计数的稳态比率，这构成了一个恒定的转录状态。然后从该比率中获取速度。
# 正速度表示基因上调，这在显示该基因未剪接 mRNA 丰度高于稳态预期水平的细胞中发生。相反，负速度表示基因下调。

# 通过设置 mode='dynamical' 可以获得完整动力学模型的解，这需要事先运行 scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='stochastic') # 计算出的速度被存储在 adata.layers 中
adata

# 然后，跨基因的速度组合可以用来估计单个细胞的未来状态。为了将速度投影到低维嵌入中，需要估计细胞间的转换概率。
# 也就是说，对于每个速度向量，我们找到与该方向一致的可能的细胞转换。转换概率是使用潜在细胞间转换和速度向量之间的余弦相关性计算的，并存储在一个表示速度图的矩阵中。得到的速度图维度为 ，
# 总结了可以通过速度向量很好地解释的可能细胞状态变化（为了提高运行速度，也可以通过设置 approx=True 在 PCA 降维空间中计算）。

scv.tl.velocity_graph(adata)
# 对于多种应用，通过应用高斯核将余弦相关性转换为实际转换概率，速度图可以转换为转换矩阵。您可以通过 scv.utils.get_transition_matrix 访问马尔可夫转换矩阵。
# 如前所述，它通过应用由 scv.tl.velocity_embedding 获取的关于转移概率的均值转移，将速度投影到低维嵌入中。
# 此外，我们可以沿着马尔可夫链追踪细胞到它们的起源和潜在命运，从而在一个轨迹中获得根细胞和端点，这是通过 scv.tl.terminal_states 获取的。

# 速度投影

# 最后，速度被投影到由 basis 指定的任何嵌入上，并以以下方式之一进行可视化
# 在细胞层面上使用 scv.pl.velocity_embedding ，
# 作为 gridlines 使用 scv.pl.velocity_embedding_grid 
# 或作为 streamlines 使用 scv.pl.velocity_embedding_stream 。

# 注意，数据已经有一个预先计算的 UMAP 嵌入和注释的簇。将应用于您自己的数据时，
# 这些可以通过 scv.tl.umap 和 scv.tl.louvain 获取。更多详情，请参阅 scanpy 教程。
# 此外，所有绘图函数默认使用 basis='umap' 和 color='clusters' ，您可以相应地设置它们。
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'clusters') # 以流线形式展示的速度矢量场可提供对发育过程的精细洞察
# 我们得到的单细胞水平速度矢量场具有最精细的分辨率，每个箭头显示单个细胞运动的方向和速度
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)


# 解释速度

# 这或许是最重要的部分，我们建议用户不要将生物学结论局限于预测速度，
# 而是通过相图(phase portraits)来检查单个基因的动态变化，以了解推断的方向如何由特定基因支持。
# 基因活性由转录调控协调。特定基因的转录诱导会导致（新转录的）未剪接前体 mRNA 增加，而相反地，抑制或缺乏转录会导致未剪接 mRNA 减少。剪接 mRNA 由未剪接 mRNA 产生，并遵循相同趋势但有时间滞后。时间是隐藏/潜变量。
# 因此，动态变化需要根据实际测量值推断：即相图中显示的剪接和未剪接 mRNA。
scv.pl.velocity(adata, ['Cpe',  'Gnao1', 'Ins2', 'Adk'], ncols=2) # 检查一些marker基因的相图

# 黑色线条对应估计的“稳态”比率，即未剪接 mRNA 与剪接 mRNA 丰度之比，该比率处于恒定的转录状态。特定基因的 RNA 速度被确定为残差，即观测值偏离该稳态线的程度。
# 正速度表明基因上调，这在未剪接 mRNA 丰度高于稳态预期值的细胞中发生。相反，负速度表明基因下调。

scv.pl.scatter(
    adata, 'Cpe', color=['clusters', 'velocity'], 
    add_outline='Ngn3 high EP, Pre-endocrine, Beta'
)
# 例如，Cpe 解释了上调的 Ngn3（黄色）到前内分泌细胞（橙色）再到β细胞（绿色）的方向性，
# 而 Adk 解释了下调的导管细胞（深绿色）到 Ngn3（黄色）再到其余内分泌细胞的方向性。

# 识别重要基因
# 我们需要一种系统的方法来识别可能有助于解释结果向量场和推断谱系的基因。
# 为此，我们可以测试哪些基因具有集群特异性差异速度表达，与剩余群体相比显著更高或更低。
# 模块 scv.tl.rank_velocity_genes 运行差异速度 t 检验，并为每个集群输出基因排名。
# 可以设置阈值（例如 min_corr ），以限制测试在选定的基因候选者上。
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

kwargs = dict(frameon=False, size=10, linewidth=1.5,
              add_outline='Ngn3 high EP, Pre-endocrine, Beta')

scv.pl.scatter(adata, df['Ngn3 high EP'][:5], ylabel='Ngn3 high EP', **kwargs)
scv.pl.scatter(adata, df['Pre-endocrine'][:5], ylabel='Pre-endocrine', **kwargs)
# 例如，Ptprs、Pclo、Pam、Abcc8、Gnas 等基因支持从 Ngn3 高表达区（黄色）到前内分泌区（橙色）再到β细胞（绿色）的方向性。


# 循环祖细胞中的速度
# RNA 速度检测到的细胞周期，通过细胞周期评分（阶段标记基因平均表达水平的标准化分数）得到生物学验证。
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])
# 对于周期性导管细胞，我们可以通过 S 期和 G2M 期标记基因进行筛选。
# 前一个模块还计算了一个 Spearman 相关系数分数，我们可以利用它对阶段标记基因进行排序/排名，然后展示它们的阶段肖像。
s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index

kwargs = dict(frameon=False, ylabel='cell cycle genes')
scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), **kwargs)

# 特别是 Hells 和 Top2a 非常适合解释循环祖细胞中的矢量场。Top2a 在 G2M 期实际达到峰值之前不久被分配了高速度。在那里，负速度完美地与紧随其后的下调相匹配。

scv.pl.velocity(adata, ['Hells', 'Top2a'], ncols=2, add_outline=True)

# 速度与一致性
# 另外两个有用的统计量：
# - 分化的速度或速率由速度向量的长度给出。
# - 向量场的相干性（即速度向量与其邻近速度的相关性）提供了置信度的度量。
scv.tl.velocity_confidence(adata)
keys = ['velocity_length', 'velocity_confidence']
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])
# 这些可以揭示细胞分化速度较慢/较快的位置，以及方向是确定/不确定的区域。
# 在聚类层面上，我们发现分化在细胞周期退出后（Ngn3 低 EP）显著加速，在β细胞产生过程中保持这一速度，而在α细胞产生过程中则减慢。


# 计算均值
df = adata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

# 速度图和伪时间
# 我们可以通过可视化速度图来描绘所有基于速度推断的细胞间连接/转换。通过设置 threshold ，可以将其限制在高概率转换上。
# 例如，该图显示了两个 Epsilon 细胞产生的阶段，分别来自早期和晚期前内分泌细胞。
scv.pl.velocity_graph(adata, threshold=.1)
# 此外，该图可用于绘制来自指定细胞的子细胞/祖细胞。在这个例子，一个前内分泌细胞被追踪到其潜在命运。

x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

# 最后，基于速度图可以计算速度伪时间。通过从图中推断根细胞分布，它测量从根细胞开始沿着图行走到达一个细胞所需的平均步数。
# 与扩散伪时间相反，它隐式地推断根细胞，并且基于有向速度图而非基于相似性的扩散核。
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')


#  PAGA 速度图
# PAGA 图抽象已被证明是轨迹推断的顶级方法。它提供了一个类似图的映射，展示了数据拓扑结构，其中加权边对应于两个簇之间的连接性。
# 在此，PAGA 通过速度推断的方向性进行了扩展。
# PAGA requires to install igraph, if not done yet.
# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='clusters') # 如果 scv.tl.paga 报错，参考https://github.com/theislab/scvelo/issues/1241，修改源代码
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
# This reads from left/row to right/column, thus e.g. assigning a confident transition from Ductal to Ngn3 low EP.
# 此表可以通过叠加到 UMAP 嵌入上的有向图进行总结
scv.pl.paga(
    adata, basis='umap', size=50, alpha=.1, 
    min_edge_width=2, node_size_scale=1.5
)