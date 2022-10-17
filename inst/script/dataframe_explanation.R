# This R scripts helps us generate data("daraframe.rda", package = "stJoincount")
# This example data has come from a Seurat object described in `inst/script/Seurat_Data_explanation.R` (the one before slimed).
# It contains the following information that is essential to this algorithm - barcode (index), cluster (they could either be categorical or numerical labels), central pixel location (imagerow and imagecol). This dataframe is simplified after combining metadata of Seurat object with spatial coordinates.

# necessary R-packages
# install.packages("Seurat")
library(Seurat)

dataPrepFromSeurat <- function(SeuratObj, label){

  metadata <- SeuratObj@meta.data
  names(metadata)[names(metadata) == label] <- "Cluster"
  coord <- GetTissueCoordinates(SeuratObj)
  merged <-  merge(metadata, coord, by = 0)
  rownames(merged) <- merged$Row.names
  df <- merged[c("imagecol", "imagerow", "Cluster")]

  return(df)
}

humanBC <- dataPrepFromSeurat(SeuratBC, "Cluster")
