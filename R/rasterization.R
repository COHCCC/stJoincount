#' Converts a labeled spatial tissue map into a raster object,
#' in which each spatial cluster is represented by a pixel coded by label assignment.
#'
#' @importFrom raster rasterize mosaic
#' @importFrom sp coordinates
#' @export
#'
#' @param sampleInfo A dataset of a human breast cancer sample containing the
#' pixel information and cluster labels for each barcode.
#'
#' @return A raster object converted from a labeled spatial tissue map.
#'
#' @examples
#' fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)

rasterizeEachCluster <- function(sampleInfo){
  sampleCoord <- sampleInfo
  sampleCoord$numericCluster <- unclass(as.factor(sampleCoord$Cluster))
  clusterNumber <- length(unique(sampleCoord$numericCluster))

  n <- extentBuffer(sampleInfo)
  r <- rasterPrep(sampleInfo, n)
  mosaicClusters <- list()

  for (i in seq_len(clusterNumber)){
    subCluster <- subset(sampleCoord, numericCluster == i)
    spdf <- subCluster
    sp::coordinates(spdf) <- c("imagerow", "imagecol") #create a Spatial object
    nam <- paste("clusterRast", i, sep = "_")
    clusterName <- assign(nam, rasterize(spdf, r, field = i, extent = jc.extent, background = 0))
    mosaicClusters <- c(mosaicClusters, clusterName)
  }
  names(mosaicClusters)[seq_len(2)] <- c('x', 'y')
  mosaicClusters$fun <- max
  mosaicClusters$na.rm <- TRUE

  mosaicIntegration <- do.call(mosaic, mosaicClusters)
  return(mosaicIntegration)
}

#' Visualization of the rasterization results and label coding of the sample.
#'
#' @importFrom raster rasterToPoints
#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_manual coord_equal theme_void theme element_blank
#' @param sampleInfo A dataset of a human breast cancer sample containing the
#' pixel information and cluster labels for each barcode.
#' @param mosaicIntegration A raster object converted from a labeled spatial tissue map.
#' @export
#'
#' @return A mosaic plot with labeled pixels.
#'
#' @examples
#' fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' mosaicIntPlot(humanBC, mosaicIntegration)
#'
mosaicIntPlot <- function(sampleInfo, mosaicIntegration){
  nameList <- customDict(sampleInfo)

  r_df <- data.frame(rasterToPoints(mosaicIntegration))
  clusterNumber <- seq_along(unique(r_df$layer))
  layerNumber <- clusterNumber - 1
  r_df$cuts <- cut(r_df$layer,breaks=layerNumber)
  r_plot <- na.omit(r_df)

  colors.raster <- c("#810505", "#f79c09", "#f30808","#f3eb17", "#bcf775", "#5cd04d", "#2b7f20", "#2cdcc2",
                     "#78cdf5", "#1b5fe4", "#811be4", "#cc6ae6", "#f59ad5", "#9c6d6d", "#0000FF")

  ggplot(data=r_plot) +
    geom_tile(aes(x=x,y=y,fill=cuts)) +
    scale_fill_manual("Cluster", values = colors.raster, labels = nameList[clusterNumber]) +
    coord_equal() +
    theme_void() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    )
}
