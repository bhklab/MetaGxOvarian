library(MetaGxOvarian)

context("Checking loadPancreasDatasets")


test_that("ensure datasets and duplicates are properly loaded from the hub and package", {
  dataAndDuplicates = MetaGxOvarian::loadOvarianDatasets()
  seData = dataAndDuplicates$summarizedExperiments
  expect_equal(is(seData[[1]])[1], "RangedSummarizedExperiment")
})
