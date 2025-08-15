library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
setwd("./scverse/")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peaks <- readPeakFile("./data/chipseq/fastq_out/Vehicle_summits.bed")
peaks1 <- list(peaks=peaks)
# promotor区间范围可以自己设定
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000) # 上下游3kb
tagMatrixList <- lapply(peaks1, getTagMatrix, windows=promoter)
#annotatePeak传入annoDb参数,可进行基因ID转换（Entrez，ENSEMBL，SYMBOL，GENENAME）
peakAnnoList <- lapply(
  peaks1, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,
  addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db"
)

tagMatrix <- getTagMatrix(peaks1$peaks, windows=promoter)

tagHeatmap(tagMatrix)

plotAvgProf(tagMatrix, xlim=c(-3000, 3000) )
?plotAvgProf

#查看peak在全基因组的位置
covplot(peaks1$peaks)
covplot(peaks1$peaks,chrs=c("chr1","chr2"))#specific chr
?covplot


peaks <- readPeakFile("./data/chipseq/fastq_out/Vehicle_summits.bed")
##用peakAnno函数对其进行peak注释，
#这里的注释主要可以分为三种，genomic注释、nearest gene注释以及peak上下游注释，
#具体解释可参照包作者写的文档peak注释。按照文章中对靶基因的要求，我们选择前两者注释即可
peakAnno <- annotatePeak(
  peaks, tssRegion=c(-3000, 3000),
  TxDb=txdb, annoDb="org.Hs.eg.db"
)
plotAnnoPie(peakAnno)
vennpie(peakAnno)

upsetplot(peakAnno)
#可视化TSS区域的TF binding loci
plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding   to TSS")

# Create a list with genes from each sample
gene = as.data.frame(peakAnno@anno@elementMetadata@listData[["SYMBOL"]])
colnames(gene)="gene"
gene1 = bitr(gene$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Hs.eg.db")
gene2 <- gene1[,2]
head(gene)
# Run GO enrichment analysis 
ego <- enrichGO(gene = gene2, 
                keyType = "ENTREZID", 
                OrgDb =org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
# Dotplot visualization
dotplot(ego, showCategory=50)
peak.n=as.data.frame(peakAnno)
write.csv()


mkdir software
cd  software/
  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToB
chmod +x bedGraphToBigWig
mkdir bin
sudo mv bedGraphToBigWig /bin
cd bin
bedGraphToBigWig
##下载fetchChromSizes，速度也很快。
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x fetchChromSizes
sudo mv fetchChromSizes  /bin
fetchChromSizes

