# 数据下载、转换、重命名
示例数据：SRA：PRJNA768891
SRR16213608
SRR16213609
SRR16213610

source /data/share/nas1/clks/config/proxy/clash-pub.proxyrc
prefetch --max-size 100G SRR16213608


# cellranger-atac

cellranger-atac 下载
wget -O cellranger-atac-2.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.2.0.tar.gz?Expires=1755222641&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=UyIonlcRf0scV~m-0YmCy0m7X9yc8avYiDdCLa8jPa6O9A~Ps-CQgGIIDqHlsGKqcGTMLzPFAS3mSbY0ScJEf2ngfAwmjECPUJIidEOiDjFEZePP~pwnnQnAvzCZr02Jd8u7hLuhJ-rIUYizlaVrWDu5wGe3kE5HQpWDTDtc-AQra8JiRSWgFVlSQY8qf5TYuoS84ufDfxh72sa4cMzJFUwLFhWrMdsG-7365IhjTXHcvcJlW~~8gXkQrGuQsQY3xMxYCA5EKVZLQ2bQHVqX7D9Ru0U7gomvD3YSRgWAZWVp~3Bh5h4SQ2rbY0sgBKJIQguffWJYc1rcprIooVF--g__"

tar -xzvf cellranger-atac-2.2.0.tar.gz
rm cellranger-atac-2.2.0.tar.gz

echo 'export PATH=/data/share/nas1/sjwlab/louhao/software/cellranger-atac-2.2.0:$PATH' >> ~/.bashrc
source ~/.bashrc 
cellranger-atac --help

# 下游

Signac：它被设计为在Seurat框架内工作，非常适合已经熟悉Seurat的用户。
这种集成允许利用Seurat强大的多模态功能，对RNA和染色质可访问性数据进行无缝联合分析。但是处理非常大规模的数据集时显得有些缓慢或消耗大量内存，特别是在数据预处理或可视化时。
#1.filtered_peak_bc_matrix.h5
#2.singlecell.csv
#3.fragments.tsv.gz
流程：创建对象 ---质控 --- LSI---降维聚类

ArchR：它是一个独立的工具，包括用于单细胞 ATAC-seq 分析的广泛功能，
例如降维、聚类和轨迹推理。它还支持与RNA-seq数据的集成，但其重点主要放在染色质可及性数据上。
#1.fragments.tsv.gz