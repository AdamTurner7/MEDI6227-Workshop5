###QCB Session 5 Workshop
# install scater as it has some nice functions to visualize and evaluate stages of your
# analyses in your workshops install.packages('remotes')
# remotes::install_github('davismcc/scater')
library(dplyr)
library(patchwork)
library(Seurat)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
dim(pbmc)
#q2 13714  2700
#(min.cells = 3 and min.features = 200)
pbmcTEST <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 2, min.features = 300)
pbmc
dim(pbmcTEST)
#q3 filters by cells with number of < 200 genes with at least 1 transcript count, #features recorded in at least these cells
#q4 view(pbmc)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#visualise data as a count matrix
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.
dense.size <- object.size(as.matrix(pbmc.data))
dense.size #data including "."
sparse.size <- object.size(pbmc.data)
sparse.size # data excluding e.g. only shows actual counts and not "."
dense.size/sparse.size
#q5 23.7 bytes

#qc data stored in metadata "percent mt"

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#once visualised set parameters to remove outliers
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#qc  should be considered by finding count peaks in number of genes, RNA depth and percent of mtgenes
#this should be considered together instead of seperately
#if distribution of covariates differ between samples, qc thresholds should be determined seperately for each sample to account for quality differences

#q6
pmbc.filtering.1 <- pbmc
pbmc.filtering.1 <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5 & percent.mt > 0.5 & nCount_RNA < 7500)
VlnPlot(pbmc.filtering.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

#q7
raw_umi_counts_loc <- pbmc@assays[["RNA"]]@layers[["counts"]]

raw_umi_counts <- as.matrix(raw_umi_counts_loc)

##investigate the effect of global normalisation by comparing the normalised data to the UMI counts data using PCA.
# export the Seurat to a scater object -
pbmc.sce <- as.SingleCellExperiment(pbmc, assay = "RNA")
detach(package:Seurat, TRUE)
# remove.packages('irlba') install.packages('irlba',type = 'source') library(irlba)
library(scater)
pbmc.sce <- scater::runPCA(pbmc.sce, exprs_values = "counts")
p1 <- scater::plotPCA(pbmc.sce, size_by = "nCount_RNA", colour_by = "nFeature_RNA") + ggtitle("PCA for UMI counts")
detach(package:scater, TRUE)


# export the Seurat to a scater object
library(Seurat)
pbmc.sce.2 <- as.SingleCellExperiment(pbmc, assay = "RNA")
detach(package:Seurat, TRUE)
library(scater)
pbmc.sce.2 <- scater::logNormCounts(pbmc.sce.2)
pbmc.sce.2 <- scater::runPCA(pbmc.sce.2, exprs_values = "logcounts")
p2 <- scater::plotPCA(pbmc.sce.2, size_by = "nCount_RNA", colour_by = "nFeature_RNA") + ggtitle("PCA for normalised data")
detach(package:scater, TRUE)

p1 + p2
