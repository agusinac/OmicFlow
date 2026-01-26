test_that("Testing omics loading from existing data objects", {
  n_cols <- 5
  n_rows <- 100
  n_vals <- n_cols * n_rows
  
  metadata <- data.table::data.table("SAMPLE_ID" = paste0("Sample_", 1:n_cols))
  features <- data.table::data.table("FEATURE_ID" = paste0("protein_", 1:n_rows))
  counts <- Matrix::Matrix(
    1:n_vals, nrow = n_rows, ncol = n_cols, 
    dimnames = list(features$FEATURE_ID, metadata$SAMPLE_ID)
  )
  
  # building up omics with metaData -> countData -> featureData
  ends_with_featureData <- omics$new(metaData = metadata)
  ends_with_featureData$countData <- counts

  expect_snapshot(ends_with_featureData)

  # building up omics with metaData -> featureData -> countData
  ends_with_countData <- omics$new(metaData = metadata)
  ends_with_countData$featureData <- features
  ends_with_countData$countData <- counts

  expect_snapshot(ends_with_countData)
})
