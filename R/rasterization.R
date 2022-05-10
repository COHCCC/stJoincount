#' Rasterize each cluster and return the integrated mosaic
#'
#' @importFrom Seurat GetTissueCoordinates
#' @import raster
#' @import spdep
#' @import sp
#' @export
#'
#' @param sample seruat object that have cluster labels attached.
#'
#' @return Integrated mosaic.
#'
#' @examples
#' fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)

rasterizeEachCluster <- function(sample){
  sampleCoord <- GetTissueCoordinates(sample)
  sampleCoord$clusters <- sample$Cluster

  n <- extentBuffer(sample)
  r <- rasterPrep(sample, n)
  mosaicClusters <- list()

  for (i in 1:length(unique(sampleCoord$clusters))){
    subCluster <- subset(sampleCoord, clusters == i)
    spdf <- subCluster
    coordinates(spdf) <- c("imagerow", "imagecol") #create a Spatial object
    nam <- paste("clusterRast", i, sep = "_")
    clusterName <- assign(nam, rasterize(spdf, r, field = i, extent = jc.extent, background = 0))
    mosaicClusters <- c(mosaicClusters, clusterName)
  }
  names(mosaicClusters)[1:2] <- c('x', 'y')
  mosaicClusters$fun <- max
  mosaicClusters$na.rm <- TRUE

  mosaicIntegration <- do.call(mosaic, mosaicClusters)
  return(mosaicIntegration)
}

#' Visulization of the rasterization results
#'
#' @import Seurat
#' @param mosaicIntegration Integrated mosaic of each cluster.
#' @export
#'
#' @return mosaic plot with integrated pixels.
#' @examples
#' fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' mosaicIntPlot(mosaicIntegration)
#'
mosaicIntPlot <- function(mosaicIntegration){
  cuts <- c(0,0.99,1.99,2.99,3.99,4.99,5.99,6.99,7.99,8.99, 9.99, 10.99, 11.99) #for setting colors per cluster
  colors.raster <- c("#FFFFFF", "#810505", "#f79c09", "#f30808","#f3eb17", "#bcf775", "#5cd04d", "#2b7f20", "#2cdcc2", "#78cdf5", "#1b5fe4", "#811be4", "#cc6ae6", "#f59ad5", "#9c6d6d")

  plot(mosaicIntegration, breaks = cuts, col = colors.raster)
  # return(p1)
}
