checkInput <- function(sampleInfo)
{
  # check overall object
  expect_type(sampleInfo, "list")

  # check colName
  cdat.cols <- c("imagerow", "imagecol", "Cluster")
  expect_true(all(cdat.cols %in% colnames(sampleInfo)))

  # check spatialCoords
  expect_true(is.numeric(sampleInfo$imagecol))
  expect_true(is.numeric(sampleInfo$imagerow))
}

test_that("the observed join counts", {
  fpath <- system.file("extdata", "dataframe.rda", package="stJoincount")
  load(fpath)
  checkInput(humanBC)

  mosaicIntegration <- rasterizeEachCluster(humanBC)
  joincount.result <- joincountAnalysis(mosaicIntegration)
  expect_type(joincount.result, "list")
})

