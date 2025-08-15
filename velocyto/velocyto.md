https://velocyto.org/velocyto.py/tutorial/cli.html#introduction

scRNAseq cellranger的输出文件,data/velocyto/results/sample1

重点关注outs文件夹，里面有我们需要的信息

drwxrwxr-x 2 louhao louhao 4.0K  8月 14 21:21 filtered_feature_bc_matrix
-rw-rw-r-- 1 louhao louhao  49M  8月 14 21:21 filtered_feature_bc_matrix.h5
-rw-rw-r-- 1 louhao louhao  669  8月 14 21:25 metrics_summary.csv
-rw-rw-r-- 1 louhao louhao 788M  8月 14 21:22 molecule_info.h5
-rw-rw-r-- 1 louhao louhao  28G  8月 14 21:16 possorted_genome_bam.bam
-rw-rw-r-- 1 louhao louhao  10M  8月 14 21:20 possorted_genome_bam.bam.bai
drwxrwxr-x 2 louhao louhao 4.0K  8月 14 21:10 raw_feature_bc_matrix
-rw-rw-r-- 1 louhao louhao  96M  8月 14 21:09 raw_feature_bc_matrix.h5
-rw-rw-r-- 1 louhao louhao 5.1M  8月 14 21:25 web_summary.html

以下是 Cell Ranger `outs` 文件夹中每个文件的详细解释：

| 文件/目录名称                     | 作用描述                                                                 | 文件类型/大小       |
|----------------------------------|--------------------------------------------------------------------------|--------------------|
| **`filtered_feature_bc_matrix`** | **过滤后的高质量细胞数据**：包含经质控过滤的细胞（真实细胞）的基因表达矩阵（barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz） | 目录 (4.0K)       |
| `filtered_feature_bc_matrix.h5`  | 上述过滤后数据的 HDF5 格式文件，适用于 Python/R 等工具直接读取                  | HDF5 文件 (49MB)   |
| `metrics_summary.csv`            | **关键指标摘要**：包含实验统计信息（如检测细胞数、中位基因数、测序饱和度等）         | CSV 文件 (669B)    |
| `molecule_info.h5`               | **分子级别原始数据**：存储每个检测到的分子（UMI）的基因组位置和细胞来源信息，用于下游深度分析 | HDF5 文件 (788MB)  |
| `possorted_genome_bam.bam`       | **比对结果文件**：测序数据比对到参考基因组的排序结果（按基因组位置排序），用于可视化或变异检测 | BAM 文件 (28GB)    |
| `possorted_genome_bam.bam.bai`   | 上述 BAM 文件的索引文件，用于快速随机访问数据                                    | BAI 索引 (10MB)    |
| **`raw_feature_bc_matrix`**      | **原始细胞数据**：包含所有检测到的条形码（包括空滴/低质量细胞）的原始基因表达矩阵      | 目录 (4.0K)       |
| `raw_feature_bc_matrix.h5`       | 上述原始数据的 HDF5 格式文件                                                  | HDF5 文件 (96MB)   |
| `web_summary.html`               | **交互式报告**：可视化实验质控结果（包括细胞数、测序深度、比对率等图表）             | HTML 报告 (5.1MB) |

### 关键说明：
1. **`filtered_` vs `raw_`**  
   - `filtered_`：仅包含通过质控的**真实细胞**（推荐用于标准分析）
   - `raw_`：包含**所有条形码**（含背景噪音，用于诊断或特殊分析）

2. **核心分析文件**  
   - `filtered_feature_bc_matrix.h5` 是下游分析（如 Seurat、Scanpy）的常用输入
   - `web_summary.html` 是快速评估数据质量的首选

3. **大文件用途**  
   - `possorted_genome_bam.bam` 用于 IGV 可视化或变异调用
   - `molecule_info.h5` 支持 Loupe 浏览器或单细胞拷贝数变异分析（如 inferCNV）

> **注**：文件大小可能因实验规模（细胞数量/测序深度）而异，其中 BAM 文件通常最大。


velocyto 分析需要
1. possorted_genome_bam.bam

2. 基因组注释文件（如果使用cellranger，那么使用cellranger配套的gtf文件）:

cp /data/share/nas1/sjwlab/louhao/reference/cellranger_ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz /data/share/nas1/sjwlab/louhao/reference/velocyto

3. 重复表达注释（gtf文件）：UCSC浏览器网页下载，然后上传到服务器

参考文件存放在reference/velocyto


# velocyto run

gzip -d /data/share/nas1/sjwlab/louhao/reference/velocyto/genes.gtf.gz

gzip -d /data/share/nas1/sjwlab/louhao/reference/velocyto/hg38_rmsk.gtf.gz

For example if we want to run the pipeline on the cellranger output folder mypath/sample01. We would do:
所以RNA速率也是每个样品进行一次分析

cd /data/share/nas1/sjwlab/louhao/scverse/data/velocyto/sample1/results/sample1/outs
samtools sort -@ 36 -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam

velocyto run10x -m /data/share/nas1/sjwlab/louhao/reference/velocyto/hg38_rmsk.gtf /data/share/nas1/sjwlab/louhao/scverse/data/velocyto/sample1/results/sample2 /data/share/nas1/sjwlab/louhao/reference/velocyto/genes.gtf >log_vel 2>&1 &

输出的数据存储在cellranger的 outs/velocyto