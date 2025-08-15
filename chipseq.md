# 下载数据
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR2242691/SRR2242691
wget 

# 将srr文件进行重命名
mv SRR8557351 IP_H3K27acVehicle
mv SRR8557352 IP_H3K27ac_XY018
mv SRR8557353 Input_H3K27ac_Vehicle
mv SRR8557354 Input_H3K27ac_XY108

# 转换为fastq文件
conda install -c bioconda parallel-fastq-dump
parallel-fastq-dump  -h
parallel-fastq-dump --sra-id Input_H3K27ac_Vehicle --threads 12  --outdir fastq_out/ --split-files --gzip
--sra-id 如果输入的是sra id，那么程序会从网上下载数据，所以我们先重命名
如果一个srr转换为两个fastq文件，说明是双端测序，一个是1，一个是2

# 数据清洗、质控（去除低质量reads、接头等）
使用 trim_galore 软件
cutadapt软件可以对NGS数据进行质量过滤，FastQC软件可以查看NGS数据的质量分布，trim_galore将两个软件封装在一起，使用起来更加的方便

conda create -n chipseq -c bioconda trim-galore
conda activate chipseq
conda install trim-galore
conda install -c bioconda bowtie2

cutadapt --version
fastqc --version


trim_galore Input_H3K27ac_Vehicle_1.fastq.gz  -q 25 --phred33 --length 25 -e 0.1 --stringency 4
trim_galore IP_H3K27acVehicle_1.fastq.gz  -q 25 --phred33 --length 25 -e 0.1 --stringency 4
trim_galore Input_H3K27ac_XY108_1.fastq.gz  -q 25 --phred33 --length 25 -e 0.1 --stringency 4
trim_galore IP_H3K27ac_XY018_1.fastq.gz  -q 25 --phred33 --length 25 -e 0.1 --stringency 4
参数具体含义看官网

## 生成质控报告
其实在生成fastq之后，先看质控报告，如果没问题就可以，如果有问题，做trim-galore，然后看质控报告
mkdir fatqc
fastqc -o fatqc/ -t 32 Input_H3K27ac_XY108_1_trimmed.fq.gz
fastqc -o fatqc/ -t 32 Input_H3K27ac_Vehicle_1_trimmed.fq.gz
fastqc -o fatqc/ -t 32 IP_H3K27ac_XY018_1_trimmed.fq.gz
fastqc -o fatqc/ -t 32 IP_H3K27acVehicle_1_trimmed.fq.gz

# 使用bowtie2进行比对
Bowtie 2 是一个超快速且内存高效的工具，用于将测序读对齐到长参考序列。它在对齐长度约为 50 到 100 个字符或 1,000 个字符的读方面表现尤为出色，尤其擅长对齐相对较长的基因组（例如哺乳动物基因组）。

## 下载参考

官网：https://bowtie-bio.sourceforge.net/bowtie2/index.shtml，提高了基因组的下载

wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip
unzip hg19.zip

参考基因组存放在：/data/share/nas1/sjwlab/louhao/index/bowtie2/ucsc.hg19

需要注意的是：genome_index指的是用于bowtie2的索引文件（如下图），而不是参考基因组本身，构建过程参考后文。genome_index需要指定路径及其共用文件名，比如我的索引文件放在/data/ref/bowtie2/hg19目录下，但是需要输入的参数为/data/ref/bowtie2/hg19/hg19。最后一个hg19指的是共用文件名。

## 正式比对
bowtie2 -h

bowtie2 -p 32 -x /data/share/nas1/sjwlab/louhao/index/bowtie2/ucsc.hg19/hg19 -U IP_H3K27acVehicle_1_trimmed.fq.gz -S IP_H3K27ac_Vehicle.sam

bowtie2 -p 32 -x /data/share/nas1/sjwlab/louhao/index/bowtie2/ucsc.hg19/hg19 -U Input_H3K27ac_Vehicle_1_trimmed.fq.gz -S Input_H3K27ac_Vehicle.sam

bowtie2 -p 32 -x /data/share/nas1/sjwlab/louhao/index/bowtie2/ucsc.hg19/hg19 -U IP_H3K27ac_XY018_1_trimmed.fq.gz -S IP_H3K27ac_XY018.sam

bowtie2 -p 32 -x /data/share/nas1/sjwlab/louhao/index/bowtie2/ucsc.hg19/hg19 -U Input_H3K27ac_XY108_1_trimmed.fq.gz -S Input_H3K27ac_XY108.sam

-x：索引文件前缀，需要指定路径和共用文件名
-p：核心数
-U：单端测序
-q：可选，输入为FASTQ文件
-f：可选，输入为FASTA文件

注意关注输出结果中的比对率，太低的话是不是测序问题，或者参考基因组选取错误

## 将sam转为bam

samtools sort -O bam -@ 32 Input_H3K27ac_Vehicle.sam >Input_H3K27ac_Vehicle.sorted.bam
samtools sort -O bam -@ 32 IP_H3K27ac_Vehicle.sam >IP_H3K27ac_Vehicle.sorted.bam
samtools sort -O bam -@ 32 IP_H3K27ac_XY018.sam >IP_H3K27ac_XY018.sorted.bam
samtools sort -O bam -@ 32 Input_H3K27ac_XY108.sam >Input_H3K27ac_XY108.sorted.bam

@：使用核心数

## 峰值的定量
macs2 callpeak -t IP_H3K27ac_Vehicle.sorted.bam -c Input_H3K27ac_Vehicle.sorted.bam -f BAM -g hs -n Vehicle -B -q 0.05
macs2 callpeak -t IP_H3K27ac_XY018.sorted.bam -c Input_H3K27ac_XY108.sorted.bam -f BAM -g hs -n XY108 -B -q 0.05

-t：IP
-c：对照
-g：基因组，人类为hs

.xls, .narrowPeak, .bed：这三个文件其实是一个意思，都是峰值位点的位置和score
.bdg文件不能直接用igv打开，可使用begGraphToBigWig、bedClip将beggraph（.bgd）转换为bigwig文件（.bw）

将IP组减去Input组，得到相减后的bdg
macs2 bdgcmp -t Vehicle_treat_pileup.bdg -c Vehicle_control_lambda.bdg -m FE -o Vehicle_FE.bdg
macs2 bdgcmp -t XY108_treat_pileup.bdg -c XY108_control_lambda.bdg -m FE -o XY108_FE.bdg

## 把bdg转换为bw文件

进入UCSC官网：https://hgdownload.soe.ucsc.edu/downloads.html

Huamn gemone--
GRCh37/hg19--Genome sequence files and select annotations (2bit, GTF, GC-content, etc)
打开后下拉到最后找到 hg19.chrom.size 文件，进行下载，如果没有txt后缀就加上去

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

bash bdg2bw.sh XY108_FE.bdg hg19.chrom.sizes.txt
bash bdg2bw.sh Vehicle_FE.bdg hg19.chrom.sizes.txt

# peaks 差异分析

在chipseq分析过程中，会使用MACS2的callpeaks命令对bam文件操作，通过比较input对照文件来找到peaks文件。但是我们的实验设计会有不同处理的样本，比如加药处理前后的样本。所以要探索一下两个不同处理的样本之间peaks有没有差异。但是callpeaks是每个处理样本与input或者igG最为control之间比较（可能有时还没有input），不能用callpeaks来做不同处理的样本之间差异分析。

目前我知道对chipseq数据进行peaks的差异分析方法有三种。
1、使用MACS2中的bdgdiff：https://github.com/macs3-project/MACS/wiki/Call-differential-binding-events
2、使用DiffBind包在R中进行
3.使用HOMER中getDifferentialPeaks命令

这里使用第三种方法



# motif分析

以MACS2软件call peak得到的bed文件分析
HOMER一个常用的motif分析软件。它通过比较两个序列集，并使用ZOOPS scoring和超几何分布（或者负二项分布）进行motif的富集分析。它主要用于ChIP-seq和promoter分析，但也可以用于核酸序列的motif分析问题。
HOMER软件可以进行多种类型的motif分析，如 promoter motif analysis，基因组位置motif分析(ChIP-seq分析中的motif分析)，利用自定义的fasta文件进行motif分析，RNA序列的motif分析（分析CLIP-seq数据中的RNAbinding elements 

HOMER进行motif分析时，需要两个数据集：
感兴趣的目标序列，如ChIP-seq分析中的peak文件
背景序列，如ChIP-seq分析中的物种全基因组序列

perl configureHomer.pl -install hg19

mkdir homer
findMotifsGenome.pl Vehicle_summits.bed hg19 homer/

Error:ls: 无法访问 '/data/share/nas1/sjwlab/louhao/conda_envs/Homer/share/homer/.//data/genomes/hg19//preparsed//hg19.*.cgbins': 没有那个文件或目录

findMotifsGeome.pl 脚本。
该脚本的输入文件支持两种类型的文件，一种是常见的bed文件格式(MACS2软件call peak得到的bed文件即可直接使用)，另一种是HOMER软件指定使用的peak文件格式。
下面详细讲一下这两种文件格式：
bed文件格式：findMotifsGeome.pl脚本用到的信息为bed文件的前六列，分别是chr,start,end,peakID,score(HOMER用不到该信息，但是必须有),strand。
peak文件格式：该文件格式需要五列信息（使用Tab分隔），分别是peakID，chr，start，end,strand。

conda create -n Homer 
conda install homer

conda info -e
cd /data/share/nas1/sjwlab/louhao/conda_envs/Homer
ls
cd share
cd homer
perl configureHomer.pl -list # 查看配置文件
perl configureHomer.pl -install human
perl configureHomer.pl -install hg19

# 找差异peak
getDifferentialPeaks -h
重要参数：
-F：差异peak富集倍数，默认是4
-P:富集的p值，默认是0.0001
-size：peak的区间长度设置

# 注释peak
既然我们已经从ChIP-seq数据中识别出了峰值，现在是时候弄清楚它们的位置和它们可能调节的基因的更多信息了。HOMER包含一个名为annotatePeaks.pl的程序，它使用peak/BED文件执行各种各样的函数。annotatePeaks.pl的用法如下：
annotatePeaks.pl -h

annotatePeaks.pl XY108_peaks.narrowPeak hg19 > scatter.txt

或者可以使用R包ChIPseeker进行注释和可视化

# 富集分析

R包 clusterprofile

# 韦恩图

intervene：一个用于交集和可视化多个基因或基因组区域集的工具

conda activate intervene

intervene --help
intervene pairwise --help
intervene venn --help
intervene upset --help

intervene venn -i Vehicle_summits.bed XY108_summits.bed -o intervene --figtype png

# 热图

deepTools: tools for exploring deep sequencing data

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.knownGene.gtf.gz

对于TSS到TES之间的区域，computeMatrix将不同基因的长度均一化到同样的长度。而对于TSS上游或者下游，则无需均一化。
computeMatrix reference-point --referencePoint TSS -b 2500 -a 2500 -S XY108_FE.bw  --regionBodyLength 5000 -R  hg19.knownGene.gtf.gz --skipZeros -o computeMatrix/XY108_FE_TSS.gz --outFileSortedRegions XY108_FE_region_genes.bed
#信号起始位点取转录前起始位点2500bp，信号终止位点取转录起始位点之后2500bp，基因body区域设置为5000bp。
##最后绘制全部peaks在TSS附近的热图和profile图
cd computeMatrix
plotHeatmap -m XY108_FE_TSS.gz -out XY108_FE_TSS.png
plotProfile --dpi 720 -m XY108_FE_TSS.gz -out XY108_FE_profile.pdf --plotFileFormat pdf

bamCoverage -p 4 -b 687.sorted.bam -o 687.sorted.bw

###如果是多个样本，可下述用computeMatrix一起读入.bw文件，再绘制图片
computeMatrix reference-point --referencePoint TSS -b 2500 -a 2500 -S ./*.bw -R ~/annotation/mm10/mm10.bed --skipZeros -o merged_TSS.gz --outFileSortedRegions merged_region_genes.bed
plotHeatmap -m merged_TSS.gz -out merged.png
plotProfile --dpi 720 -m merged_TSS.gz -out merged_profile.pdf --plotFileFormat pdf

