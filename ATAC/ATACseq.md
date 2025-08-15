# 数据下载

source /data/share/nas1/clks/config/proxy/clash-pub.proxyrc
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15175276/SRR15175276

# 数据转换
conda activate sra-env

mv SRR15175276 atac
mkdir fastq_out
parallel-fastq-dump --sra-id atac --threads 16  --outdir fastq_out/ --split-files --gzip
--sra-id 如果输入的是sra id，那么程序会从网上下载数据，所以我们先重命名
如果一个srr转换为两个fastq文件，说明是双端测序，一个是1，一个是2

cd fastq_out

# 数据清洗、质控（去除低质量reads、接头等）
使用 trim_galore 软件
cutadapt软件可以对NGS数据进行质量过滤，FastQC软件可以查看NGS数据的质量分布，trim_galore将两个软件封装在一起，使用起来更加的方便

conda activate chipseq

trim_galore atac_1.fastq.gz  -q 25 --phred33 --length 25 -e 0.1 --stringency 4

参数具体含义看官网

## 生成质控报告
其实在生成fastq之后，可以先看质控报告，如果没问题就可以，如果有问题，做trim-galore，然后再看质控报告
mkdir fastqc
fastqc -o fastqc/ -t 32 atac_1.fastq.gz    # trim galore之前
fastqc -o fastqc/ -t 32 atac_1_trimmed.fq.gz

可用multiqc综合所有qc
使用multiqc来把所有的质量评估放在一起观察
multqc -o <out.path> *.fastqc.zip

# 比对

mkdir hg19_genome
cd hg19_genome
wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
tar -zxvf hg19_genome.tar.gz

其实atacseq更常用bowtie2，但是这里为了演示，使用了hisat2，STAR是转为RNAseq数据分析设计的，而单细胞RNAseq有cellranger一站式分析，空间转录组有spaceranger
hisat2 -p 16 -x /data/share/nas1/sjwlab/louhao/scverse/data/ATACseq/fastq_out/hg19_genome/hg19/genome -U atac_1_trimmed.fq.gz -S atac.sam


hisat官网下载index（不自己构建的话）
参考
index文件直接去官网下载：http://daehwankimlab.github.io/hisat2/download/


hisat用户手册：http://ccb.jhu.edu/software/hisat2/manual.shtml
更适合RNA-seq，因为索引策略的使用，可以相对轻松比对跨区域的read（可变剪切）
Tophat2耗时久，STAR吃内存。hisat兼有两个的耗时、内存消耗的优点
索引的构建
注意，如果要用到gtf中exon和split site信息，要先用hisat2自带的py软件生成，然后添加到参数--exon 和参数--ss后
hisat2-build [options]<reference_in> <ht2_index_base>hisat2 [options]* -x <ht2-idx> {-1 <m1> -2<m2>|-U<r>}[-S
#参数：<sam>]
#-P线程数
#<reference_in>参考基因组fasta文件list，如果为list，使用逗号分开
#<ht2_base>索引文件的前缀名
#--ss>与--exon一起使用，提供拼接位点信息。只支持hisat2自己格式
#（需要用hisat2_extract_splice_sites.py与 gtf来生成）
#--exon指定外显子列表（只支持hisat2自己格式，需要通过
hisat2_extract_exons.py与gtf文件来生成）从GTF提取剪接位点。

# sam 格式转换

sam格式转换为bam：
samtools view -Sb atac.sam>atac.bam
samtools view -h
BAM是SAM文件的二进制版本。它存储与SAM文件相同的比对数据，但以二进制格式存储，从而节省存储空间并提高处理效率。

排序：
samtools sort -@ 16 atac.bam -o atac.sort.bam

索引bam文件 
samtools index atac.sort.bam #生成文件默认是原文件名+.bai也可以指定在原文件后面直接指定

BAM 文件的索引文件（通常为 .bai 文件）是用于 快速访问 该文件中的数据。例如，许多可视化工具（如 IGV）依赖于索引文件来快速加载和显示比对数据。如果没有索引，加载数据时会非常缓慢，甚至无法加载。

# sambamba 去除重复（可选）

conda activate sambamba

sambamba markdup -r atac.bam atac.rmdup.bam

samtools view -h -f2 -q 30 atac.rmdup.bam |grep -v chrM| samtools sort -O bam -@ 16 -o - > atac.last.bam
bedtools bamtobed -i atac.last.bam > atac.bed

samtools和bedtools在同一个conda环境中
# 峰值定量
mkdir ../peaks
macs2 callpeak -t atac.bed -g mm --nomodel --shift -100 --extsize 200 -n atac --outdir ../peaks/ # 这句话报错，用下面的，所以好像没必要转换为bed文件

macs2 callpeak -t atac.sort.bam -g mm --nomodel --shift -100 --extsize 200 -n atac --outdir ../peaks/

# bam To bw
conda activate deeptools
bamCoverage -b atac.sort.bam -o atac.bw

# peak 注释

同chipseq