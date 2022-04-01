#' Finding optimal number of buffer for extent.
#'
#' @import Seurat
#' @import raster
#' @import spdep
#' @import sp
#' @export
#'
#' @param sample seruat object that have cluster labels attached.
#'
#' @return optimal number of buffer for extent
#'
#' @examples
#' data('humanBC')
#' n <- extentBuffer(humanBC)

extentBuffer <- function(sample){
  sampleCoord <- GetTissueCoordinates(sample)
  sampleCoord$clusters <- sample$Cluster
  clusterNumber <- length(unique(sampleCoord$clusters))

  n = 10
  for (j in 1:20){
    k = 0
    r <- rasterPrep(sample, n)
    for (i in 1:clusterNumber){
      subCluster <- subset(sampleCoord, clusters == i)
      spdf <- subCluster
      coordinates(spdf) <- c("imagerow", "imagecol") #create a Spatial object
      nam <- paste("clusterRast", i, sep = "_")
      clusterName <- assign(nam, rasterize(spdf, r, field = i, extent = jc.extent, background = 0))
      valueCheck <- sum(clusterName@data@values)/i
      if (valueCheck != nrow(subCluster)){
        # print("The number of pixels is not equal to the number of barcodes")
        break
      } else { k = k + 1}
    }
    if (k == clusterNumber){
      print(paste("Found optimal number of buffer. n =", n, sep = " "))
      return(n)
      break
    }
    n = n + 5
  }
  print("No optimal number of buffer found, using n = 110 instead.")
  print("In this case, there may be minor deviations in the subsequent calculation process. The results are for reference only")
  return(n)
}
