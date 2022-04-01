#' Finding the extent with buffer applied.
#'
#' @import Seurat
#' @import raster
#' @export
#'
#' @param sample seruat object that have cluster labels attached.
#' @param n buffer for extent
#'
#' @return raster layer with calculated resolution and extent with buffer applied
#'
#' @examples
#' data('humanBC')
#' raster <- rasterPrep(humanBC, 10)

rasterPrep <- function(sample, n){
  #get coordinates
  sampleCoord <- GetTissueCoordinates(sample)
  sampleCoord$clusters <- sample$Cluster
  coordSummary <- as.data.frame(apply(sampleCoord,2,summary))

  #create raster that will apply to all subsets
  imagerow.min <- as.integer(coordSummary$imagerow[1])
  imagerow.max <- as.integer(coordSummary$imagerow[6])
  imagecol.min <- as.integer(coordSummary$imagecol[1])
  imagecol.max <- as.integer(coordSummary$imagecol[6])
  jc.extent <- extent(imagerow.min-n, imagerow.max+n, imagecol.min-n, imagecol.max+n)

  resolutionList <- resolutionCalc(sample)
  r <- raster(resolution = resolutionList, ext = jc.extent)
  return(r)
}
