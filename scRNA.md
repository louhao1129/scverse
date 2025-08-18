# 数据下载、转换、重命名
tmux new -s sra
prefetch --max-size 100G SRR13177101 > log_01 2>&1 &
prefetch --max-size 100G SRR13177102 > log_02 2>&1 &
prefetch --max-size 100G SRR13177103 > log_03 2>&1 &
prefetch --max-size 100G SRR13177104 > log_04 2>&1 &



转换为fastq
parallel-fastq-dump --sra-id atac --threads 32  --outdir fastq_out/ --split-files --gzip

parallel-fastq-dump --sra-id SRR13177101.sra --threads 32  --outdir sample1/ --split-files --gzip
parallel-fastq-dump --sra-id SRR13177102.sra --threads 32  --outdir sample1/ --split-files --gzip
parallel-fastq-dump --sra-id SRR13177103.sra --threads 32  --outdir sample1/ --split-files --gzip
parallel-fastq-dump --sra-id SRR13177104.sra --threads 32  --outdir sample1/ --split-files --gzip

cd sample1
ls -lh  # 查看文件大小，根据大小，进行重命名sample1是自己起的，后面的内容是sra数据库中作者提供的，后面的是给软件读的

mv SRR13177101_1.fastq.gz sample1_S13_L001_I1_001.fastq.gz
mv SRR13177101_2.fastq.gz sample1_S13_L001_R1_001.fastq.gz
mv SRR13177101_3.fastq.gz sample1_S13_L001_R2_001.fastq.gz

mv SRR13177102_1.fastq.gz sample1_S13_L002_I1_001.fastq.gz
mv SRR13177102_2.fastq.gz sample1_S13_L002_R1_001.fastq.gz
mv SRR13177102_3.fastq.gz sample1_S13_L002_R2_001.fastq.gz

mv SRR13177103_1.fastq.gz sample1_S13_L003_I1_001.fastq.gz
mv SRR13177103_2.fastq.gz sample1_S13_L003_R1_001.fastq.gz
mv SRR13177103_3.fastq.gz sample1_S13_L003_R2_001.fastq.gz

mv SRR13177104_1.fastq.gz sample1_S13_L004_I1_001.fastq.gz
mv SRR13177104_2.fastq.gz sample1_S13_L004_R1_001.fastq.gz
mv SRR13177104_3.fastq.gz sample1_S13_L004_R2_001.fastq.gz

# cellranger

cellranger教程：https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials

下载cellranger：https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links
tar -xzvf cellranger-9.0.1.tar.gz
rm cellranger-9.0.1.tar.gz

echo 'export PATH=/data/share/nas1/sjwlab/louhao/software/cellranger-9.0.1:$PATH' >> ~/.bashrc
source ~/.bashrc 

下载参考基因组
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

注：我还下载了cellranger8，可以使用
/data/share/nas1/sjwlab/louhao/software/cellranger-8.0.1/cellranger -h


mkdir results
cd results

cellranger count --help # 默认是使用的v9

cellranger count --id=sample1 \
                   --transcriptome=/data/share/nas1/sjwlab/louhao/reference/cellranger_ref/refdata-gex-GRCh38-2024-A \
                   --fastqs=/data/share/nas1/sjwlab/louhao/scverse/data/velocyto/sample1 \
                   --sample=sample1 \
                   --create-bam=true \
                   --localcores=32 \
                   --localmem=200 \
                   --nosecondary



Outputs:
- Run summary HTML:                         /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/web_summary.html
- Run summary CSV:                          /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/metrics_summary.csv
- BAM:                                      /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/possorted_genome_bam.bam
- BAM BAI index:                            /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/possorted_genome_bam.bam.bai
- BAM CSI index:                            null
- Filtered feature-barcode matrices MEX:    /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            null
- Per-molecule read information:            /data/share/nas1/sjwlab/louhao/scverse/data/velocity/sample1/results/sample1/outs/molecule_info.h5
- CRISPR-specific analysis:                 null
- Antibody aggregate barcodes:              null
- Loupe Browser file:                       null
- Feature Reference:                        null
- Target Panel File:                        null
- Probe Set File:                           null
- cell_types:                               null
- Symlink web_summary_cell_types.html:      null

Alerts:
Not running cell annotation as secondary analysis has been disabled!

Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

2025-08-14 21:25:30 Shutting down.
Saving pipestance info to "sample1/sample1.mri.tgz"

# 同法处理剩下样本

source /data/share/nas1/clks/config/proxy/clash-pub.proxyrc

prefetch --max-size 100G SRR13177105 > log_05 2>&1 &
prefetch --max-size 100G SRR13177106 > log_06 2>&1 &
prefetch --max-size 100G SRR13177107 > log_07 2>&1 &
prefetch --max-size 100G SRR13177108 > log_08 2>&1 &

转换为fastq
parallel-fastq-dump --sra-id atac --threads 32  --outdir fastq_out/ --split-files --gzip

mkdir -p sample2
parallel-fastq-dump --sra-id SRR13177105.sra --threads 32  --outdir sample2/ --split-files --gzip
parallel-fastq-dump --sra-id SRR13177106.sra --threads 32  --outdir sample2/ --split-files --gzip
parallel-fastq-dump --sra-id SRR13177107.sra --threads 32  --outdir sample2/ --split-files --gzip
parallel-fastq-dump --sra-id SRR13177108.sra --threads 32  --outdir sample2/ --split-files --gzip

cd sample2
ls -lh  # 查看文件大小，根据大小，进行重命名sample1是自己起的，后面的内容是sra数据库中作者提供的，后面的是给软件读的

mv SRR13177105_1.fastq.gz sample2_S12_L001_I1_001.fastq.gz
mv SRR13177105_2.fastq.gz sample2_S12_L001_R1_001.fastq.gz
mv SRR13177105_3.fastq.gz sample2_S12_L001_R2_001.fastq.gz

mv SRR13177106_1.fastq.gz sample2_S12_L002_I1_001.fastq.gz
mv SRR13177106_2.fastq.gz sample2_S12_L002_R1_001.fastq.gz
mv SRR13177106_3.fastq.gz sample2_S12_L002_R2_001.fastq.gz

mv SRR13177107_1.fastq.gz sample2_S12_L003_I1_001.fastq.gz
mv SRR13177107_2.fastq.gz sample2_S12_L003_R1_001.fastq.gz
mv SRR13177107_3.fastq.gz sample2_S12_L003_R2_001.fastq.gz

mv SRR13177108_1.fastq.gz sample2_S12_L004_I1_001.fastq.gz
mv SRR13177108_2.fastq.gz sample2_S12_L004_R1_001.fastq.gz
mv SRR13177108_3.fastq.gz sample2_S12_L004_R2_001.fastq.gz


mkdir results
cd results

cellranger count --help # 默认是使用的v9

cellranger count --id=sample2 \
                   --transcriptome=/data/share/nas1/sjwlab/louhao/reference/cellranger_ref/refdata-gex-GRCh38-2024-A \
                   --fastqs=/data/share/nas1/sjwlab/louhao/scverse/data/velocyto/sample2 \
                   --sample=sample2 \
                   --create-bam=true \
                   --localcores=32 \
                   --localmem=200 \
                   --nosecondary