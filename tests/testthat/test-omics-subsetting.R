## Load example data
test <- proteomics$new(
  metaData = "input/proteomics/metadata.csv",
  countData = "input/proteomics/counts.csv"
)

test_that("`omics$sample_subset()` -- Argument checks", {
  expect_error(test$sample_subset(gender > 5), "Failed to evaluate subset expression `gender > 5`: object 'gender' not found")
  expect_error(test$sample_subset(test1), "Failed to evaluate subset expression `test1`: object 'test1' not found")

  expect_snapshot(test$sample_subset(CONTRAST_treatment == "nothing"), error = TRUE)
})

test_that("`omics$sample_subset()` -- Behavioral checks", { 
  # Perform metadata subset
  test$reset()
  test$sample_subset(CONTRAST_treatment == "tumor")
  expect_snapshot(test)

  # Test active binding
  test$reset()
  expect_snapshot(test)
  expect_snapshot(test$metaData)
  expect_snapshot(test$countData)
  expect_snapshot(test$treeData)

  test$reset()
  test$countData <- test$countData[1:50, ]
  expect_snapshot(test)
})

test_that("`omics$feature_subset()` -- Argument checks", {
  expect_error(test$feature_subset(gender > 5), "Failed to evaluate subset expression `gender > 5`: object 'gender' not found")
  expect_error(test$feature_subset(test1), "Failed to evaluate subset expression `test1`: object 'test1' not found")

  expect_snapshot(test$feature_subset(grepl("nothing", FEATURE_ID)), error = TRUE)
})

test_that("`omics$feature_subset()` -- Behavioral checks", { 
  # Pass a grepl to subset
  test$reset()
  test$feature_subset(grepl("3|4", FEATURE_ID))
  expect_snapshot(test)

  # active binding subset
  test$reset()
  test$featureData
  test$featureData <- test$featureData[1:200, ]
  expect_snapshot(test)
})

test_that("`omics$samplepair_subset()` -- Argument checks", {
  test$metaData[, "SAMPLEPAIR_ID" := NULL]
  expect_snapshot(test$samplepair_subset(), error = TRUE)
  expect_snapshot(test$samplepair_subset(num_unique_pairs = 2), error = TRUE)
})

test_that("`omics$samplepair_subset()` -- Behavioral checks", { 
  test$reset()
  test$samplepair_subset()
  expect_snapshot(test)
  expect_equal(unique(test$metaData$SAMPLEPAIR_ID), c("S001", "S002", "S004"))
})