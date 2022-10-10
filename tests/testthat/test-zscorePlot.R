test_that("Visualization of the z-score matrix", {
  fpath <- system.file("extdata", "dataframe.rda", package="stJoincount")
  load(fpath)

  mosaicIntegration <- rasterizeEachCluster(humanBC)
  joincount.result <- joincountAnalysis(mosaicIntegration)
  matrix <- zscoreMatrix(humanBC, joincount.result)
  p <- zscorePlot(matrix)
  expect_s3_class(p, "pheatmap")
})
