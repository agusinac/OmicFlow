test_that("Testing `feature_merge`", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  features_original <- data.table::copy(taxa$featureData)
  counts_original <- taxa$countData
  
  ## Feature merge
  taxa$feature_merge(
    feature_rank = "Genus",
    feature_filter = c("uncultured")
  )

  expect_equal(
    Matrix::colSums(counts_original[features_original[Genus == "Woesearchaeales"]$FEATURE_ID, ]),
    taxa$countData[taxa$featureData[Genus == "Woesearchaeales"]$FEATURE_ID, ] 
  )
})