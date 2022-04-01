#' Automatically calculating the resolution for rasterization.
#'
#' @importFrom stats median
#' @import graphics
#' @import Seurat
#' @export
#'
#' @param sample seruat object that have cluster labels attached.
#'
#' @return length and height of resolution.
#'
#' @examples
#' data('humanBC')
#' resolutionList <- resolutionCalc(humanBC)

resolutionCalc <- function(sample){
  copyCoord <- GetTissueCoordinates(sample)
  copyCoord$int.row <- as.integer(copyCoord$imagerow)
  copyCoord$int.col <- as.integer(copyCoord$imagecol)

  diffCol <- list()
  diffRow <- list()

  for (i in 1:nrow(copyCoord)){
    a <- copyCoord$int.row[i]
    targetCol <- subset(copyCoord, int.row == a)[order(subset(copyCoord, int.row == a)$imagecol),]
    subtractCol <- diff(as.matrix(targetCol$imagecol))
    outliersA <- boxplot(subtractCol, plot=FALSE)$out
    if (length(outliersA) > 0){
      subtractCol <- subtractCol[-which(subtractCol %in% outliersA),]
      differenceInCol <- mean(subtractCol)/2
    } else {
      differenceInCol <- mean(subtractCol)/2
    }
    diffCol <- c(diffCol, differenceInCol)

    b <- copyCoord$int.col[i]
    targetRow <- subset(copyCoord, int.col == b)[order(subset(copyCoord, int.col == b)$imagerow),]
    subtractRow <- diff(as.matrix(targetRow$imagerow))
    outliersB <- boxplot(subtractRow, plot=FALSE)$out
    if (length(outliersB) > 0){
      subtractRow <- subtractRow[-which(subtractRow %in% outliersB),]
      differenceInRow <- mean(subtractRow)
    } else {
      differenceInRow <- mean(subtractRow)
    }
    diffRow <- c(diffRow, differenceInRow )
  }

  X <- lapply(diffRow, function(x) x[!is.na(x)])
  X <- median(unlist(X))
  Y <- lapply(diffCol, function(x) x[!is.na(x)])
  Y <- median(unlist(Y))

  # checking the orientation
  if (abs(X-Y) < 2){
    X = X/2
    Y = Y*2
  }

  return(c(X, Y))
}
