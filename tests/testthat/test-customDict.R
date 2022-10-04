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

test_that("create label dictionary", {
  fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
  load(fpath)
  checkInput(humanBC)

  dict <- customDict(humanBC)
  expect_type(dict, "list")
})
