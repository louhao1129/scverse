# 数据下载、重命名、转换

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR9137271/SRR9137271

mv SRR9137271 hcc

conda activate sra-env
parallel-fastq-dump -h

mkdir fastq_out
parallel-fastq-dump --sra-id hcc --threads 16  --outdir fastq_out/ --split-files --gzip
--sra-id 如果输入的是sra id，那么程序会从网上下载数据，所以我们先重命名
如果一个srr转换为两个fastq文件，说明是双端测序，一个是1，一个是2



# 过滤reads：trim_galore
可以先看质控报告
Trim Galore是对FastQC和Cutadapt的包装。适用于所有高通过滤标准量测序，包括RRBS(Reduced Representation Bisulfite-Seq）,Illumina、Nextera和smallRNA测序平台的双端和单端数据。主要功能包括两步：
第一步首先去除低质量碱基，然后去除3'末端的adapter,如果没有指定具体的adapter，程序会自动检测前1million的序trim_galore列，然后对比前12-13bp的序列是否符合以下类型的adapter:

过滤标准
1、一般的过滤标准：
2、去除含有接头的reads
3、过滤去除低质量值的reads
4、去除含有N比例大于5%的read

conda activate chipseq
trim_galore hcc_1.fastq.gz  -q 25 --phred33 --length 25 -e 0.1 --stringency 4


# QC

其实在生成fastq之后，先看质控报告，如果没问题就可以，如果有问题，做trim-galore，然后看质控报告
mkdir fatqc
fastqc -o fatqc/ -t 32 Input_H3K27ac_XY108_1_trimmed.fq.gz

软件
fastqc可以对*.bam*.sam *.fq*.fq.gz进行质量评估。
fastqc可以通过-t指定多线程操作，多线程是同时处理多个输入文件，几个线程可以同时处理几个文件，
单个文件使用多线程似乎没有意义

multiqc综合所有qc
使用multiqc来把所有的质量评估放在一起观察
multqc -o <out.path> *.fastqc.zip

# 比对：bowtie2、hisat2、STAR

把reads回贴到基因组上面

https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes
    Parent Directory                               -   
    hg19.ensGene.gtf.gz       2020-01-10 09:45   26M  # ensemble
    hg19.knownGene.gtf.gz     2024-04-25 10:23   35M  # symbol
    hg19.ncbiRefSeq.gtf.gz    2021-05-17 10:35   19M  # ncbi
    hg19.refGene.gtf.gz       2020-01-10 09:45   21M  # ref


wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.knownGene.gtf.gz


使用bowtie2进行比对：
bowtie2 -p 32 -x /data/share/nas1/sjwlab/louhao/index/bowtie2/ucsc.hg19/hg19 -U hcc_1_trimmed.fq.gz -S hcc.sam




chipseq我们使用了bowtie2，使用的是uscs官网上的索引，这里我们使用hisat2，并且学习如何下载基因组并且创建index
进入ucsc官网：https://genome.ucsc.edu/，找到download--Genome Data，进入，
我们选择下载人类hg19作为示例数据：Genome sequence files and select annotations (2bit, GTF, GC-content, etc)，patient Directory，bigZip
下载基因组文件：wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
解压缩：tar -zxvf chromFa.tar.gz
解压出来的新文件都是同一类型：chr*.fa。
chr1.fa
chr10.fa
ochr11.fa
chr11_gl000202_random.fa
chr12.fa
chr13.fa
[..]
cat *.fa>hg19.fa #将解压出来的所有.fa结尾的文件合并为一个文件，即hg19.fa，到此为止，参考基因组就准备好了。
建立HISAT2索引文件：hisat2-build -p 35 hg19.fa genome


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

# sam格式转换

sam格式转换samtools
sam格式：The Sequence Alignment／Map format，存储了序列比对情况的文件，是一个文本文件，都以tab分列。体积较大。分为头部区和主体区。
官方文档：http://samtools.github.io/hts-specs/SAMvl.pdf

排序和转BAM
samtools官方文档：http://www.htslib.org/doc/samtools.html
samtools默认使用坐标排序
当未指定samtool的索引类型时，默认是使用BAI索引。
详细使用方法可以使用man samtoolssort查看

sam格式转换为bam：samtools view -Sb hcc.sam>hcc.bam
samtools view -h
BAM是SAM文件的二进制版本。它存储与SAM文件相同的比对数据，但以二进制格式存储，从而节省存储空间并提高处理效率。

samtools sort -@ 16 hcc.bam -o hcc.sort.bam

索引bam文件 
samtools index sample.sort.bam #生成文件默认是原文件名+.bai也可以指定在原文件后面直接指定

# 定量(计数)

已经得到了某基因回贴到的reads数，要转换为count数目

由于序列读长的限制和基因组的同源性，一条reads可能比对到多个基因上，而且基因之间也存在overlap,在对这些特殊情况进行处理时，htseq-count内置了以下3种模式
union
intersection-strict
intersection-nonempty
通过--mode参数指定某种模式，默认值为union。
这3种模式在判断一条reads是否属于某个feature时，有不同的判别标准

在运行速度上，featurecounts比htseq-count快很多倍，而且feature-count不仅支持基因/转录本的定量，也支持exon等单个feature的定量，所以更推荐使用featurecounts来定量。
和featurecounts一样，htseq-count也是一款进行raw count定量的软件。该软件采用python语言进行开发，集成在HTseq这个包中。
HTSeq提供了许多处理NGS数据的功能，htseq-count只是其中进行定量分析的一个模块。
HTSeq使用注意事项
1、HTSeq是对有参考基因组的转录组测序数据进行表达量分析的，其输入文件必须有SAM和GTF文件。
2、一般情况下HTSeq得到的Counts结果会用于下一步不同样品间的基因表达量差异分析，而不是一个样品内部基因的表达量比较。因此，HTSeq设置了-a参数的默认值10，来忽略掉比对到多个位置的reads信息，其结果有利于后续的差异分析。
3、输入的GTF文件中不能包含可变剪接信息，否则HTSeq会认为每个可变剪接都是单独的基因，导致能比对到多个可变剪接转录本上的reads的计算结果是ambiguous，从而不能计算到基因的count中。即使设置-i参数的值为transcript _id，其结果一样是不准确的，只是得到transcripts的表达量。
htseq-count的设计思想和featurecounts非常类似，也包含了feature和meta-features两个概念。
对于转录组数据而言，feature指的是exon,而meta-feature可以是gene,也可以是transcript。
进行定量分析需要以下两个文件
1、比对的BAM/SAM文件
2、基因组的GTF文件
对于双端数据，要求输入sort之后的BAM文件。



htseq-count 命令参数
-f|--format default:sam 设置输入文件的格式，该参数的值可以是sam或bam。
-r|--order default:name 设置sam或bam文件的排序方式，该参数的值可以是name或pos。前者表示按read名进行排序，后者表示按比对的参考基因组位置进行排序。若测序数据是双末端测序，当输入sam/bam文件是按pos方式排序的时候，两端reads的比对结果在sam/bam文件中一般不是紧邻的两行，程序会将reads对的第一个比对结果放入内存，直到读取到另一端read的比对结果。因此，选择pos可能会导致程序使用较多的内存，它也适合于未排序的sam/bam文件。而pos排序则表示程序认为双末端测序的reads比对结果在紧邻的两行上，也适合于单端测序的比对结果。很多其它表达量分析软件要求输入的sam/bam文件是按pos排序的，但HTSeq推荐使用name排序，且一般比对软件的默认输出结果也是按name进行排序的。
-s|--stranded default:yes 设置是否是链特异性测序。该参数的值可以是yes,no或reverse。no表示非链特异性测序；若是单端测序，yes表示read比对到了基因的正义链上；若是双末端测序，yes表示readl比对到了基因正义链上，read2比对到基因负义链上；reverse表示双末端测序情况下与yes值相反的结果。根据说明文件的理解，一般情况下
双末端链特异性测序，该参数的值应该选择reverse（本人暂时没有测试该参数）。
htseq-count -f bam -s no -r name -i gene_id -t exon hcc.sort.bam /data/share/nas1/sjwlab/louhao/scverse/data/RNAseq/fastq_out/hg19.knownGene.gtf.gz> hcc.count.txt

