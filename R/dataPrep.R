#' Data Preparation from Seurat Object
#'
#' This data preparation function creates a data.frame from a Seurat Object.
#' It extracts the pixel information and cluster labels of each barcode from user's input Seurat Object
#' and generate a data.frame with a certain format which is required for the algorithm.
#' If the user has customized labels, this function will change the column name to "Cluster" when generating
#' the data.frame to make it consistent to the required format.
#'
#' @importFrom Seurat GetTissueCoordinates
#' @export
#'
#' @param SeuratObj input Seurat object that contains labels for each barcode
#' @param label the column name of the label information in "meta.data"
#'
#' @return A data.frame contains the pixel information and cluster labels for each barcode the sample.
#' The index contains barcodes, and at least three other columns that have these information are required and the column names should be the same as following:
#' "imagerow": The row pixel coordinate of the center of the spot
#' "imagecol": The column pixel coordinate of the center of the spot
#' "Cluster": The label that corresponding to this barcode
#'
#' @examples
#' fpath <- system.file("extdata", "SeuratBC.rda", package="stJoincount")
#' load(fpath)
#' df <- dataPrepFromSeurat(seuratBC, "Cluster")


dataPrepFromSeurat <- function(SeuratObj, label){

  metadata <- SeuratObj@meta.data
  names(metadata)[names(metadata) == label] <- "Cluster"
  coord <- GetTissueCoordinates(SeuratObj)
  merged <-  merge(metadata, coord, by = 0)
  rownames(merged) <- merged$Row.names
  df <- merged[c("imagecol", "imagerow", "Cluster")]

  return(df)
}

#' Data Preparation from SpatialExperiment Object
#'
#' This data preparation function creates a data.frame form Spatial Experiment Object.
#' It extracts the pixel information and cluster labels of each barcode from user's input Seurat Object
#' and generate a data.frame with a certain format which is required for the algorithm.
#' If the user has customized labels, this function will change the column name to "Cluster" when generating
#' the data.frame to make it consist to the required format.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @export
#'
#' @param SpeObj input SpatialExperiment object that contains labels for each barcode
#' @param label the column name of the label information in "colData"
#'
#' @return A data.frame contains the pixel information and cluster labels for each barcode of the sample.
#' The index contains barcodes, and at least three other columns that have these information are required and the column names should be the same as following:
#' "imagerow": The row pixel coordinate of the center of the spot
#' "imagecol": The column pixel coordinate of the center of the spot
#' "Cluster": The label that corresponding to this barcode
#'
#' @examples
#' fpath <- system.file("extdata", "SpeBC.rda", package="stJoincount")
#' load(fpath)
#' df <- dataPrepFromSpE(SpeObjBC, "label")


dataPrepFromSpE <- function(SpeObj, label){

  metadata <- colData(SpeObj)
  coord <- spatialCoords(SpeObj)
  merged <-  merge(metadata, coord, by = 0)
  rownames(merged) <- merged$Row.names

  merged <- merged %>% rename("imagecol" = "pxl_col_in_fullres",
                              "imagerow" = "pxl_row_in_fullres",
                              "Cluster" = label)

  df <- merged[c("imagecol", "imagerow", "Cluster")]

  return(df)
}
