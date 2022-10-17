# This R scripts helps us generate data("SeuratBC", package = "stJoincount")
# User needs to directory containing the "filtered_feature_bc_matrix.h5"(HDF5), and three subdirectory "spatial"(Sptial imaging data) , "analysis"(Clustering analysis), and "filtered_feature_bc_matrix"(barcode matrix(filtered)).
# These data can be found:https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0

# necessary R-packages
# install.packages("Seurat")
library(Seurat)

spatialDataPrep <- function(filename){
  inputSample <- Load10X_Spatial(filename)
  sampleMeta <- read.csv(paste(filename, "/analysis/clustering/graphclust/clusters.csv", sep = ""), sep = ',')
  if (colnames(sampleMeta)[1] == "barcode"){
    colnames(sampleMeta)[colnames(sampleMeta) == "barcode"] <- "Barcode"
    colnames(sampleMeta)[colnames(sampleMeta) == "cluster"] <- "Cluster"
  }
  rownames(sampleMeta) <- sampleMeta$Barcode
  sampleCluster <- AddMetaData(object = inputSample, metadata = sampleMeta)
  Idents(sampleCluster) <- "Cluster"
  return(sampleCluster)
}

# folder pathway that contains feature-bc-matrix, analysis and spatial information as previous described
humanBC <- spatialDataPrep("~/human_breast_cancer")

# Following steps used for sliming data to meet the size requirement of bioconductor packages
slim.BC <- DietSeurat(
  humanBC,
  counts = FALSE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL
)
idx <- sample(ncol(slim.BC),20)
seuratBC <- slim.BC[, idx]

