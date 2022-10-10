test_that("generate data.frame from Seurat object", {
  fpath <- system.file("extdata", "SeuratBC.rda", package="stJoincount")
  load(fpath)

  df <- dataPrepFromSeurat(SeuratBC)
  expect_type(df, "data.frame")
})

test_that("generate data.frame from Spatial Experiment object", {
  fpath <- system.file("extdata", "SpeBC.rda", package="stJoincount")
  load(fpath)

  df <- dataPrepFromSpE(SpeObjBC)
  expect_type(df, "data.frame")
})
