#' When sample coordinates finds a suitable buffer to ensure that each cluster is not lost
#' in the process of converting the spot to pixel,
#' apply this buffer to this function to find a suitable rectangle for the rasterlayer
#'
#' @importFrom raster raster extent
#' @export
#'
#' @param sampleInfo A dataset of a human breast cancer sample containing the
#' pixel information and cluster labels for each barcode.
#' @param n buffer for extent (from function extentBuffer).
#'
#' @return raster layer with calculated resolution and extent with buffer applied
#'
#' @examples
#' fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' raster <- rasterPrep(humanBC, 15)

rasterPrep <- function(sampleInfo, n){
  imagerow.min <- as.integer(min(sampleInfo$imagerow))
  imagerow.max <- as.integer(max(sampleInfo$imagerow))
  imagecol.min <- as.integer(min(sampleInfo$imagecol))
  imagecol.max <- as.integer(max(sampleInfo$imagecol))
  jc.extent <- extent(imagerow.min-n, imagerow.max+n, imagecol.min-n, imagecol.max+n)

  resolutionList <- resolutionCalc(sampleInfo)
  r <- raster(resolution = resolutionList, ext = jc.extent)
  return(r)
}

