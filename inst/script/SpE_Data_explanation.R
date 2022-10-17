# This R scripts helps us generate data("SpeBC", package = "stJoincount")
# These data can be found:https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0
# These steps are already explained in the tutorial at:
# https://lmweber.org/OSTA-book/feature-selection.html

# necessary R-packages
# install.packages("SpatialExperiment")
# install.packages("STexampleData")

library(SpatialExperiment)
library(STexampleData)

spe <- Visium_humanDLPFC()

# QUALITY CONTROL (QC)
library(scater)

# subset to keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
# calculate per-spot QC metrics
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# select QC thresholds
qc_lib_size <- colData(spe)$sum < 600
qc_detected <- colData(spe)$detected < 400
qc_mito <- colData(spe)$subsets_mito_percent > 28
qc_cell_count <- colData(spe)$cell_count > 10
# combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
colData(spe)$discard <- discard
# filter low-quality spots
spe <- spe[, !colData(spe)$discard]

# NORMALIZATION

library(scran)
# calculate logcounts using library size factors
spe <- logNormCounts(spe)

# FEATURE SELECTION

# remove mitochondrial genes
spe <- spe[!is_mito, ]
# fit mean-variance relationship
dec <- modelGeneVar(spe)
# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)

# DIMENSIONALITY REDUCTION

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)
# compute UMAP on top 50 PCs
spe <- runUMAP(spe, dimred = "PCA")
# update column names
#colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
# graph-based clustering
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)

# Following steps used for sliming data to meet the size requirement of bioconductor packages
#remove both images since we do not need them in this package
imgData(spe <- rmvImg(spe))
imgData(spe <- rmvImg(spe))
set.seed(123)
idx <- sample(ncol(spe), 20)
SpeObjBC <- spe[, idx]

