test_that("Testing omics loading", {
  metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
  counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
  features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

  with_features <- omics$new(
    countData = counts_file,
    featureData = features_file,
    metaData = metadata_file
  )
  expect_snapshot(with_features)
  
  without_features <- omics$new(
    countData = counts_file,
    metaData = metadata_file
  )
  expect_snapshot(without_features) 
})
