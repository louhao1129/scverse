#BiocManager::install("pcaMethods")
#library(devtools)
#install_github("velocyto-team/velocyto.R")
#remotes::install_github("satijalab/seurat-wrappers")
library(pcaMethods)
library(velocyto.R)
library(SeuratWrappers)
library(Seurat)
getwd()
setwd("./scverse/")

# curl::curl_download(
#   url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', 
#   destfile = './data/velocity/SCG71.loom'
# )
# Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE")
# Sys.setenv(RHDF5_USE_FILE_LOCKING="FALSE")
ldat <- ReadVelocity(file = "./data/velocity/SCG71.loom")
bm <- as.Seurat(x = ldat)
bm <- SCTransform(object = bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)


bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
