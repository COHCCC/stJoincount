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

test_that("cxreate the raster object converted from a labeled spatial tissue map", {
  fpath <- system.file("script", "humanBC.rda", package="stJoincount")
  load(fpath)
  checkInput(humanBC)

  expectedRaster <- rasterizeEachCluster(humanBC)
  expect_s4_class(expectedRaster[[1]], "RasterLayer")
})

