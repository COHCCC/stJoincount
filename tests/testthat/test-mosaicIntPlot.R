test_that("Visualization of the rasterization results", {
  fpath <- system.file("script", "humanBC.rda", package="stJoincount")
  load(fpath)

  mosaicIntegration <- rasterizeEachCluster(humanBC)
  p <- mosaicIntPlot(humanBC, mosaicIntegration)
  expect_s3_class(p, "ggplot")
})
