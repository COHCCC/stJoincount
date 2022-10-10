#' When we create the rasterlayer, there will be a rectangular range.
#' It is often necessary to provide a buffer to ensure that subsequent functions
#' do not result in blank or missed pixels.
#' This function is to find the right buffer for the sample coordinates
#' so that each cluster is not lost in the process of converting a spot to a pixel.
#'
#' @importFrom raster raster
#' @importFrom sp coordinates
#' @export
#'
#' @param sampleInfo A data.frame contains the pixel information and cluster labels for each barcode of a human breast cancer sample.
#' The index contains barcodes, and at least three other columns that have these information are required and the column names should be the same as following:
#' "imagerow": The row pixel coordinate of the center of the spot
#' "imagecol": The column pixel coordinate of the center of the spot
#' "Cluster": The label that corresponding to this barcode
#'
#' @return optimal number of buffer for extent
#'
#' @examples
#' fpath <- system.file("extdata", "dataframe.rda", package="stJoincount")
#' load(fpath)
#' n <- extentBuffer(humanBC)

extentBuffer <- function(sampleInfo){
  sampleCoord <- sampleInfo
  sampleCoord$numericCluster <- unclass(as.factor(sampleCoord$Cluster))
  clusterNumber <- length(unique(sampleCoord$numericCluster))

  n <- 10
  for (j in seq_len(20)){
    k <- 0
    r <- rasterPrep(sampleInfo, n)
    for (i in seq_len(clusterNumber)){
      subCluster <- subset(sampleCoord, numericCluster == i)
      spdf <- subCluster
      sp::coordinates(spdf) <- c("imagerow", "imagecol")
      nam <- paste("clusterRast", i, sep = "_")
      clusterName <- assign(nam, rasterize(spdf, r, field = i, extent = jc.extent, background = 0))
      valueCheck <- sum(clusterName@data@values)/i
      if (valueCheck != nrow(subCluster)){
        break
      } else { k <- k + 1}
    }
    if (k == clusterNumber){
      message(paste("Find optimal number of n. n =", n, sep = " "))
      return(n)
      break
    }
    n <- n + 5
  }
  message("No optimal number found, using n = 110 instead.")
  message("In this case, there may be minor deviations in the subsequent calculation process.
        The results are for reference only")
  return(n)
}

