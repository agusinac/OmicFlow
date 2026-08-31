test_that("`omics` -- Argument checks", {

  ## Ensuring `metaData` is supplied
  expect_snapshot(omics$new(), error = TRUE)
  expect_snapshot(omics$new(featureData = features_file), error = TRUE)
  expect_snapshot(omics$new(countData = counts_sparse_file), error = TRUE)
  expect_snapshot(omics$new(metaData = data.frame()), error = TRUE)
  expect_snapshot(omics$new(metaData = data.table::data.table()), error = TRUE)

  ## Checking errors
  expect_snapshot(omics$new(metaData = metadata_file, featureData = data.frame()), error = TRUE)
  expect_snapshot(omics$new(metaData = metadata_file, featureData = data.table::data.table()), error = TRUE)

  expect_snapshot(omics$new(metaData = metadata_file, countData = data.frame()), error = TRUE)
  expect_snapshot(omics$new(metaData = metadata_file, countData = data.table::data.table()), error = TRUE)
  expect_snapshot(omics$new(metaData = metadata_file, countData = matrix(0)), error = TRUE)
})

test_that("`omics` -- Behavioral checks", { 
  # Loading only metaData
  expect_snapshot(
    omics$new(
      metaData = metadata_file
    )
  )

  # Loading metaData and featureData
  expect_snapshot(
    omics$new(
      metaData = metadata_file,
      featureData = features_file
    )
  )

  # Loading metaData and countData without rownames
  expect_snapshot(
    test <- omics$new(
      metaData = metadata_file,
      countData = counts_sparse_file
    )
  )
  expect_equal(rownames(test$countData)[1], "feature_1")
  expect_equal(rownames(test$countData), test$featureData$FEATURE_ID)
  expect_equal(colnames(test$countData), test$metaData$SAMPLE_ID)

  # Loading metaData and countData with rownames
  expect_snapshot(
    test <- omics$new(
      metaData = metadata_file,
      countData = counts_sparse_with_rownames_file
    )
  )
  expect_equal(rownames(test$countData)[1], "GTGTCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAGACAGACGTGAAATCCCCGGGCTTAACCTGGGAACTGCGTTTGTGACTGCAAGGCTAGAGTACGGCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGATATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGTCGATACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGAAACCCTAGTAGTCC")
  expect_equal(rownames(test$countData), test$featureData$FEATURE_ID)
  expect_equal(colnames(test$countData), test$metaData$SAMPLE_ID)

  # Intialising `countData` from file (mixed with features)
  metadata <- data.table::data.table(SAMPLE_ID = c("S1", "S2", "S3"))
  counts <- data.table::data.table(S1 = c(2,3,0), S2 = c(2,"",1), S3 = c(2,1,NA), proteins = c("ZEB1", "ZEB2", "MAPK"))

  # Intialising `countData` from data.frame (mixed with features)
  expect_no_error(
    test <- omics$new(
      metaData = metadata,
      countData = as.data.frame(counts)
    )
  )
  expect_snapshot(test)
  expect_equal(rownames(test$countData), test$featureData$FEATURE_ID)
  expect_equal(test$featureData$proteins, counts$proteins)
  expect_equal(test$countData[3,1], 0) # Real zero
  expect_equal(test$countData[2,2], 0) # String replaced by zero
  expect_equal(test$countData[3,3], 0) # NA replaced by zero
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

  # Intialising `countData` from data.table (mixed with features)
  expect_no_error(
    test <- omics$new(
      metaData = metadata,
      countData = counts
    )
  )
  expect_snapshot(test)
  expect_equal(rownames(test$countData), test$featureData$FEATURE_ID)
  expect_equal(test$featureData$proteins, counts$proteins)
  expect_equal(test$countData[3,1], 0) # Real zero
  expect_equal(test$countData[2,2], 0) # String replaced by zero
  expect_equal(test$countData[3,3], 0) # NA replaced by zero
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

  # From filepath
  tmp_metadata_file <- paste0(tempdir(), "tmp_metadata.csv")
  tmp_count_file <- paste0(tempdir(), "tmp_counts.csv")

  data.table::fwrite(x = metadata, file = tmp_metadata_file)
  data.table::fwrite(x = counts, file = tmp_count_file)

  expect_no_error(
    test <- omics$new(
      metaData = tmp_metadata_file,
      countData = tmp_count_file
    )
  )
  expect_snapshot(test)
  expect_equal(rownames(test$countData), test$featureData$FEATURE_ID)
  expect_equal(test$featureData$proteins, counts$proteins)
  expect_equal(test$countData[3,1], 0) # Real zero
  expect_equal(test$countData[2,2], 0) # String replaced by zero
  expect_equal(test$countData[3,3], 0) # NA replaced by zero
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")
  ## remove files
  file.remove(tmp_metadata_file)
  file.remove(tmp_count_file)

  # Loading all three components with sparse counts
  expect_snapshot(
    sparse <- omics$new(
      metaData = metadata_file,
      featureData = features_file,
      countData = counts_sparse_file
    )
  )
  expect_equal(rownames(sparse$countData), sparse$featureData$FEATURE_ID)
  expect_s4_class(sparse$countData, "sparseMatrix")
  expect_s3_class(sparse$metaData, "data.table")
  expect_s3_class(sparse$featureData, "data.table")

  # Loading all three components with dense counts
  expect_snapshot(
    dense <- omics$new(
      metaData = metadata_file,
      featureData = features_file,
      countData = counts_dense_file
    )
  )
  expect_equal(rownames(dense$countData), dense$featureData$FEATURE_ID)
  expect_s4_class(dense$countData, "sparseMatrix")
  expect_s3_class(dense$metaData, "data.table")
  expect_s3_class(dense$featureData, "data.table")

  ## Checking active bindings and sync behavior
  #-------------------------------------------------------------
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
  expect_equal(rownames(ends_with_featureData$countData), ends_with_featureData$featureData$FEATURE_ID)
  expect_equal(colnames(ends_with_featureData$countData), ends_with_featureData$metaData$SAMPLE_ID)
  expect_s4_class(ends_with_featureData$countData, "sparseMatrix")
  expect_s3_class(ends_with_featureData$metaData, "data.table")
  expect_s3_class(ends_with_featureData$featureData, "data.table")
  expect_snapshot(ends_with_featureData)

  # building up omics with metaData -> featureData -> countData
  ends_with_countData <- omics$new(metaData = metadata)
  ends_with_countData$featureData <- features
  ends_with_countData$countData <- counts
  expect_equal(rownames(ends_with_countData$countData), ends_with_countData$featureData$FEATURE_ID)
  expect_equal(colnames(ends_with_countData$countData), ends_with_countData$metaData$SAMPLE_ID)
  expect_s4_class(ends_with_countData$countData, "sparseMatrix")
  expect_s3_class(ends_with_countData$metaData, "data.table")
  expect_s3_class(ends_with_countData$featureData, "data.table")
  expect_snapshot(ends_with_countData)
})
