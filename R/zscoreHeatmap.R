#' Z-score matrix.
#'
#' @import Seurat
#' @import reshape2
#' @export
#'
#' @param sample seruat object that have cluster labels attached
#' @param joincount.result calculated result from join count analysis
#'
#' @return z-score matrix
#'
#' @examples
#' data('humanBC')
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' joincount.result <- joincountAnalysis(mosaicIntegration)
#' matrix <- zscoreMatrix(humanBC, joincount.result)

zscoreMatrix <- function(sample, joincount.result){
  clusterNumbers <- length(unique(sample$Cluster))
  jcMatrix <- data.frame(matrix(NA, nrow = clusterNumbers, ncol = clusterNumbers))
  colnames(jcMatrix) <- as.character(c(1:clusterNumbers))
  rownames(jcMatrix) <- as.character(c(1:clusterNumbers))

  for (i in 1:clusterNumbers){
    for (j in 1:clusterNumbers){
      index1 <- paste(i, ':', j, sep = "")
      jcMatrix[i,j] <- joincount.result[index1, 'z-value']
      jcMatrix[j,i] <- joincount.result[index1, 'z-value']
    }
  }
  jcMatrix <- round(jcMatrix, 2)
  return(jcMatrix)
}


#' Z-score heatmap.
#'
#' @import Seurat
#' @import pheatmap
#' @export
#'
#' @param zscoreMatrix calculated and reshaped z-score matirx from join count analysis
#'
#' @return Heatmap plot
#' @examples
#' data('humanBC')
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' joincount.result <- joincountAnalysis(mosaicIntegration)
#' matrix <- zscoreMatrix(humanBC, joincount.result)
#' zscorePlot(matrix)

zscorePlot <- function(zscoreMatrix){
  pheatmap(zscoreMatrix, legend = TRUE, name = "z-score",
                   cluster_rows=FALSE, cluster_cols=FALSE,
                   display_numbers = TRUE, fontsize_number = 9,
                   column_names_side = c("top"), row_names_side = c("left"), angle_col = c("0"),
                   show_colnames = TRUE, show_rownames = TRUE, number_color = "black")
}
