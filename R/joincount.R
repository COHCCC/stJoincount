#' Join count analysis.
#'
#' @import Seurat
#' @import spdep
#' @importFrom stats na.omit
#' @export
#'
#' @param mosaicIntegration RasterLayer that output from Function rasterization.
#'
#' @return Join count analysis results.
#'
#' @examples
#' fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' joincount.result <- joincountAnalysis(mosaicIntegration)

joincountAnalysis <- function(mosaicIntegration){
  nbList <- cell2nb(nrow = nrow(mosaicIntegration), ncol = ncol(mosaicIntegration), type = "queen")

  emptyPos <- which(mosaicIntegration@data@values == 0)
  subList <- subset.nb(nbList,
                       !(1:length(nbList) %in% emptyPos))
  temp.subList <- lapply(subList, as.character)
  noNeighborPos <- which(sapply(temp.subList, FUN =  function(X) "0" %in% X))

  filteredList <- subset.nb(subList,
                            !(1:length(subList) %in% noNeighborPos))

  weightsList <- nb2listw(filteredList,
                          style = "B",
                          zero.policy = NULL)

  fx <- (as.factor(mosaicIntegration@data@values))[-emptyPos]
  if (length(noNeighborPos) > 0){
    fx <- fx[-noNeighborPos]
  }

  if (isTRUE(length(fx) == length(weightsList$neighbours) && length(fx) == length(weightsList$weights))){
    joincount.filtered <- joincount.multi(fx = fx,
                                          listw = weightsList,
                                          zero.policy = FALSE,
                                          adjust.n = FALSE)
  } else {
    print("The length of mosaiced data and the length of weights list is not equal.")
  }
  joincount.result <- as.data.frame(joincount.filtered)
  joincount.result <- na.omit(joincount.result)

  return(joincount.result)
}

