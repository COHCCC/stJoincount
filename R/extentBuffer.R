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
#' @param sampleInfo A dataset of a human breast cancer sample containing the
#' pixel information and cluster labels for each barcode
#'
#' @return optimal number of buffer for extent
#'
#' @examples
#' fpath <- system.file("script", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' n <- extentBuffer(humanBC)

extentBuffer <- function(sampleInfo){
  sampleCoord <- sampleInfo
  sampleCoord$numericCluster <- unclass(as.factor(sampleCoord$Cluster))
  clusterNumber <- length(unique(sampleCoord$numericCluster))

  n = 10
  for (j in seq_len(20)){
    k = 0
    r <- rasterPrep(sampleInfo, n)
    for (i in seq_len(clusterNumber)){
      subCluster <- subset(sampleCoord, numericCluster == i)
      spdf <- subCluster
      sp::coordinates(spdf) <- c("imagerow", "imagecol") #create a Spatial object
      nam <- paste("clusterRast", i, sep = "_")
      clusterName <- assign(nam, rasterize(spdf, r, field = i, extent = jc.extent, background = 0))
      valueCheck <- sum(clusterName@data@values)/i
      if (valueCheck != nrow(subCluster)){
        # print("The number of pixels is not equal to the number of barcodes")
        break
      } else { k = k + 1}
    }
    if (k == clusterNumber){
      message(paste("Find optimal number of n. n =", n, sep = " "))
      return(n)
      break
    }
    n = n + 5
  }
  message("No optimal number found, using n = 110 instead.")
  message("In this case, there may be minor deviations in the subsequent calculation process.
        The results are for reference only")
  return(n)
}

