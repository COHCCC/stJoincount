#' Prepare input data (seurat object), add label to the object
#'
#' @import Seurat
#' @importFrom utils read.csv
#' @export
#'
#' @param filename folder pathway that contains feature-bc-mstrix, analysis and spatial information.
#'
#' @return seruat  object that have cluster labels attached.
#'
spatialDataPrep <- function(filename){
  inputSample <- Load10X_Spatial(filename)
  sampleTransform <- SCTransform(inputSample, assay = "Spatial", verbose = FALSE)
  sampleMeta <- read.csv(paste(filename, "/analysis/clustering/graphclust/clusters.csv", sep = ""), sep = ',')
  if (colnames(sampleMeta)[1] == "barcode"){
    colnames(sampleMeta)[colnames(sampleMeta) == "barcode"] <- "Barcode"
    colnames(sampleMeta)[colnames(sampleMeta) == "cluster"] <- "Cluster"
  }
  rownames(sampleMeta) <- sampleMeta$Barcode
  sampleCluster <- AddMetaData(object = sampleTransform, metadata = sampleMeta)
  Idents(sampleCluster) <- "Cluster"
  return(sampleCluster)
}


