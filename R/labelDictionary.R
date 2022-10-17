#' Generate dict of cluster names
#'
#' Create a dictionary with categorical cluster labels as values
#' and their converted numerical labels as keys
#'
#' @param sampleInfo A data.frame contains the pixel information and cluster labels for each barcode of a human breast cancer sample.
#' The index contains barcodes, and at least three other columns that have these information are required and the column names should be the same as following:
#' "imagerow": The row pixel coordinate of the center of the spot
#' "imagecol": The column pixel coordinate of the center of the spot
#' "Cluster": The label that corresponding to this barcode
#' @export
#'
#' @return A dictionary with categorical cluster labels as values
#' and their converted numerical labels as keys.
#'
#' @examples
#' fpath <- system.file("extdata", "dataframe.rda", package="stJoincount")
#' load(fpath)
#' nameList <- customDict(humanBC)
#'
customDict <- function(sampleInfo){
  sampleCoord <- sampleInfo
  sampleCoord$numericCluster <- unclass(as.factor(sampleCoord$Cluster))
  clusterNumber <- length(unique(sampleCoord$numericCluster))

  nameList <- list()
  key <- sampleCoord$numericCluster
  value <- sampleCoord$Cluster
  for (i in seq_len(length(key))){
    nameList[key[i]] <- value[i]
  }

  return(nameList)
}
