## Load example data
metadata_file <- "input/proteomics/metadata.csv"
counts_with_rownames_file <- "input/proteomics/counts.csv"
counts_without_rownames_file <- "input/proteomics/counts_without_rownames.csv"
tree_file <- "input/proteomics/tree.newick"

test_that("`proteomics` -- Argument checks", {

  ## Ensuring `metaData` is supplied
  expect_snapshot(proteomics$new(), error = TRUE)
  expect_snapshot(proteomics$new(countData = counts_with_rownames_file), error = TRUE)
  expect_snapshot(proteomics$new(metaData = data.frame()), error = TRUE)
  expect_snapshot(proteomics$new(metaData = data.table::data.table()), error = TRUE)

  ## Checking errors
  expect_snapshot(proteomics$new(metaData = metadata_file, featureData = data.frame()), error = TRUE)
  expect_snapshot(proteomics$new(metaData = metadata_file, featureData = data.table::data.table()), error = TRUE)

  expect_snapshot(proteomics$new(metaData = metadata_file, countData = data.frame()), error = TRUE)
  expect_snapshot(proteomics$new(metaData = metadata_file, countData = data.table::data.table()), error = TRUE)
  expect_snapshot(proteomics$new(metaData = metadata_file, countData = matrix(0)), error = TRUE)

  expect_snapshot(proteomics$new(metaData = metadata_file, treeData = data.frame()), error = TRUE)
  expect_snapshot(proteomics$new(metaData = metadata_file, treeData = ape::rtree(50)), error = TRUE)
})

test_that("`proteomics` -- Behavioral checks", { 
  # Loading counts file with rownames
  expect_snapshot(
    test <- proteomics$new(
      metaData = metadata_file,
      countData = counts_with_rownames_file
    )
  )
  expect_equal(all(rownames(test$countData) == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(colnames(test$countData) == test$metaData$SAMPLE_ID), TRUE)
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

  # Loading counts file without rownames
  expect_snapshot(
    test <- proteomics$new(
      metaData = metadata_file,
      countData = counts_without_rownames_file
      )
  )
  expect_equal(all(rownames(test$countData) == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(colnames(test$countData) == test$metaData$SAMPLE_ID), TRUE)
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

   # Loading with tree
  expect_snapshot(
    test <- proteomics$new(
      metaData = metadata_file,
      countData = counts_with_rownames_file,
      treeData = tree_file
    )
  )
  expect_equal(all(rownames(test$countData) == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(rownames(test$countData) == test$treeData$tip.label), TRUE)
  expect_equal(all(test$treeData$tip.label == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(colnames(test$countData) == test$metaData$SAMPLE_ID), TRUE)
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

  # Checking loading proteomics from pre-loaded test
  expect_snapshot(
    prot_ref <- proteomics$new(
      countData = test$countData,
      metaData = test$metaData,
      treeData = test$treeData,
      featureData = test$featureData
    )
  )
  expect_equal(all(rownames(test$countData) == rownames(prot_ref$countData)), TRUE)
  expect_equal(all(colnames(test$countData) == colnames(prot_ref$countData)), TRUE)
  expect_equal(all(test$featureData$FEATURE_ID == prot_ref$featureData$FEATURE_ID), TRUE)
  expect_equal(all(test$treeData$tip.label == prot_ref$treeData$tip.label), TRUE)
  expect_equal(all(test$metaData$SAMPLE_ID == prot_ref$metaData$SAMPLE_ID), TRUE)
  expect_equal(inherits(test$countData, "sparseMatrix"), inherits(prot_ref$countData, "sparseMatrix"))
  expect_equal(class(test$metaData)[1], class(prot_ref$metaData)[1])
  expect_equal(class(test$featureData)[1], class(prot_ref$featureData)[1])
})