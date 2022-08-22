#' This function provides a heatmap of z-scores resulting
#' from the join count analysis for all possible label pairs.
#'
#' @param sampleInfo A dataset of a human breast cancer sample containing the
#' pixel information and cluster labels for each barcode.
#' @param joincount.result calculated result from join count analysis
#' @export
#'
#' @return z-score matrix resulting from the join count analysis for all possible label pairs
#'
#' @examples
#' fpath <- system.file("script", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' joincount.result <- joincountAnalysis(mosaicIntegration)
#' matrix <- zscoreMatrix(humanBC, joincount.result)

zscoreMatrix <- function(sampleInfo, joincount.result){
    clusterNumbers <- length(unique(sampleInfo$Cluster))
    jcMatrix <- data.frame(matrix(NA, nrow = clusterNumbers, ncol = clusterNumbers))
    nameList <- customDict(sampleInfo)
    clusterName <- nameList[seq_len(clusterNumbers)]
    colnames(jcMatrix) <- clusterName
    rownames(jcMatrix) <- clusterName

    for (i in seq_len(clusterNumbers)){
      for (j in seq_len(clusterNumbers)){
        index1 <- paste(i, ':', j, sep = "")
        jcMatrix[i,j] <- joincount.result[index1, 'z-value']
        jcMatrix[j,i] <- joincount.result[index1, 'z-value']
      }
    }
    jcMatrix <- round(jcMatrix, 2)
    return(jcMatrix)
}

#' Visulization of Z-score heatmap.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @export
#'
#' @param zscoreMatrix calculated and reshaped z-score matirx from join count analysis.
#'
#' @return Heatmap plot
#'
#' @examples
#' fpath <- system.file("script", "humanBC.rda", package="stJoincount")
#' load(fpath)
#' mosaicIntegration <- rasterizeEachCluster(humanBC)
#' joincount.result <- joincountAnalysis(mosaicIntegration)
#' matrix <- zscoreMatrix(humanBC, joincount.result)
#' zscorePlot(matrix)

zscorePlot <- function(zscoreMatrix){
    pheatmap(zscoreMatrix, legend = TRUE, name = "z-score",
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE, fontsize_number = 9,
             color=colorRampPalette(c("navy", "white", "yellow", "orange", "red"))(50),
             column_names_side = c("top"), row_names_side = c("left"), angle_col = c("0"),
             show_colnames = TRUE, show_rownames = TRUE, number_color = "black")
}

