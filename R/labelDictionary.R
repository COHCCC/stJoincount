#' Create a dictionary with categorical cluster labels as values
#' and their converted numerical labels as keys
#'
#' @param sampleInfo A dataset of a human breast cancer sample containing the
#' pixel information and cluster labels for each barcode.
#' @export
#'
#' @return A mosaic plot with labeled pixels.
#' @examples
#' fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
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
