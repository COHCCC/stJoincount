#' Join count analysis
#'
#' This function performes multi-categorical join count analysis of
#' the rasterized sample. A neighbors list is then created from the rasterized sample,
#' in which adjacent and diagonal neighbors for each pixel are identified.
#'
#' @importFrom spdep cell2nb nb2listw joincount.multi subset.nb
#' @importFrom stats na.omit
#' @export
#'
#' @param mosaicIntegration A raster object converted from a labeled spatial tissue map from Function rasterization.
#'
#' @return A data.frame that contains the observed join counts, the expected count under conditions of spatial randomness, the variance calculated under non-free sampling, and calculated Z-score.
#'
#' @examples
#' fpath <- system.file("extdata", "dataframe.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' joincount.result <- joincountAnalysis(mosaicIntegration)

joincountAnalysis <- function(mosaicIntegration){
  nbList <- cell2nb(nrow = nrow(mosaicIntegration), ncol = ncol(mosaicIntegration), type = "queen")

  emptyPos <- which(mosaicIntegration@data@values == 0)
  subList <- subset.nb(nbList,
                       !(seq_len(length(nbList)) %in% emptyPos))
  temp.subList <- lapply(subList, as.character)
  noNeighborPos <- which(vapply(temp.subList, FUN =  function(X) "0" %in% X, logical(1)))

  filteredList <- subset.nb(subList,
                            !(seq_len(length(subList)) %in% noNeighborPos))

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
    stop("The length of mosaiced data and the length of weights list is not equal.")
  }
  joincount.result <- as.data.frame(joincount.filtered)
  joincount.result <- na.omit(joincount.result)

  return(joincount.result)
}
