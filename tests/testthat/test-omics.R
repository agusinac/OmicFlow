test_that("Testing omics loading from file", {
  metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
  counts_sparse_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
  counts_sparse_with_rownames_file <- system.file("extdata", "counts_with_rownames.tsv", package = "OmicFlow")
  counts_dense_file <- system.file("extdata", "counts_dense.tsv", package = "OmicFlow")
  features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

  # Loading data with features and with sparse or dense formats from file
  with_features_sparse <- omics$new(
    countData = counts_sparse_file,
    featureData = features_file,
    metaData = metadata_file
  )
  expect_snapshot(with_features_sparse)

  with_features_dense <- omics$new(
    countData = counts_dense_file,
    featureData = features_file,
    metaData = metadata_file
  )
  expect_snapshot(with_features_dense)

  # Loading data without features but with/without rownames from counts file
  without_features_with_rownames <- omics$new(
    countData = counts_sparse_with_rownames_file,
    metaData = metadata_file
  )
  expect_snapshot(without_features_with_rownames)
  
  without_features_without_rownames <- omics$new(
    countData = counts_file,
    metaData = metadata_file
  )
  expect_snapshot(without_features_without_rownames)
})
