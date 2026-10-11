## Load example data
biom_hdf5 <- "input/metagenomics/biom_with_taxonomy_hdf5.biom"
biom_json <- "input/metagenomics/biom_with_taxonomy_json.biom"
tree <- "input/metagenomics/rooted_tree.newick"

test_that("`metagenomics` -- Argument checks", {

  ## Ensuring `metaData` is supplied
  expect_snapshot(metagenomics$new(), error = TRUE)
  expect_snapshot(metagenomics$new(biomData = "nonexisting.biom"), error = TRUE)
  expect_snapshot(metagenomics$new(featureData = features_file), error = TRUE)
  expect_snapshot(metagenomics$new(countData = counts_sparse_file), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = data.frame()), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = data.table::data.table()), error = TRUE)

  ## Checking errors
  expect_snapshot(metagenomics$new(metaData = metadata_file, featureData = data.frame()), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = metadata_file, featureData = data.table::data.table()), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = metadata_file, biomData = "nonexisting.biom"), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = metadata_file, biomData = metadata_file), error = TRUE)

  expect_snapshot(metagenomics$new(metaData = metadata_file, countData = data.frame()), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = metadata_file, countData = data.table::data.table()), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = metadata_file, countData = matrix(0)), error = TRUE)

  expect_snapshot(metagenomics$new(metaData = metadata_file, treeData = data.frame()), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = metadata_file, treeData = ape::rtree(50)), error = TRUE)
  expect_snapshot(metagenomics$new(metaData = metadata_file, biomData = biom_hdf5, treeData = ape::rtree(50)), error = TRUE)
})

test_that("`metagenomics` -- Behavioral checks", { 
  # Loading biom hdf5
  expect_snapshot(
    test <- metagenomics$new(
      metaData = metadata_file,
      biomData = biom_hdf5
    )
  )
  expect_equal(all(rownames(test$countData) == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(colnames(test$countData) == test$metaData$SAMPLE_ID), TRUE)
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

  # Loading biom json
  expect_snapshot(
    test <- metagenomics$new(
      biomData = biom_json,
      metaData = data.table::data.table(SAMPLE_ID = c(
        "Sample1", "Sample2", "Sample3",
        "Sample4", "Sample5", "Sample6"))
      )
  )
  expect_equal(all(rownames(test$countData) == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(colnames(test$countData) == test$metaData$SAMPLE_ID), TRUE)
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

   # Loading biom hdf5 with tree
  expect_snapshot(
    test <- metagenomics$new(
      metaData = metadata_file,
      biomData = biom_hdf5,
      treeData = tree
    )
  )
  expect_equal(all(rownames(test$countData) == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(rownames(test$countData) == test$treeData$tip.label), TRUE)
  expect_equal(all(test$treeData$tip.label == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(colnames(test$countData) == test$metaData$SAMPLE_ID), TRUE)
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

  # Checking loading metagenomics from pre-loaded test
  expect_snapshot(
    taxa_ref <- metagenomics$new(
      countData = test$countData,
      metaData = test$metaData,
      treeData = test$treeData,
      featureData = test$featureData
    )
  )
  expect_equal(all(rownames(test$countData) == rownames(taxa_ref$countData)), TRUE)
  expect_equal(all(colnames(test$countData) == colnames(taxa_ref$countData)), TRUE)
  expect_equal(all(test$featureData$FEATURE_ID == taxa_ref$featureData$FEATURE_ID), TRUE)
  expect_equal(all(test$treeData$tip.label == taxa_ref$treeData$tip.label), TRUE)
  expect_equal(all(test$metaData$SAMPLE_ID == taxa_ref$metaData$SAMPLE_ID), TRUE)
  expect_equal(inherits(test$countData, "sparseMatrix"), inherits(taxa_ref$countData, "sparseMatrix"))
  expect_equal(class(test$metaData)[1], class(taxa_ref$metaData)[1])
  expect_equal(class(test$featureData)[1], class(taxa_ref$featureData)[1])

  ## Handles `feature_names` with extra labels that are not present
  expect_snapshot(
    test <- metagenomics$new(
      metaData = metadata_file,
      biomData = biom_hdf5,
      feature_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "variants")
    )
  )
  expect_equal(all(rownames(test$countData) == test$featureData$FEATURE_ID), TRUE)
  expect_equal(all(colnames(test$countData) == test$metaData$SAMPLE_ID), TRUE)
  expect_s4_class(test$countData, "sparseMatrix")
  expect_s3_class(test$metaData, "data.table")
  expect_s3_class(test$featureData, "data.table")

  ## Gives warning if less `feature_names` are provided than found in the `featureData`
  expect_snapshot(
    test <- metagenomics$new(
      metaData = metadata_file,
      biomData = biom_hdf5,
      feature_names = c("Kingdom", "Phylum", "Class")
    )
  )
})