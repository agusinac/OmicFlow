test_that("Testing Metagenomics reading and writing of BIOM files", {
  # Loading hdf5 format
  taxa_hdf5 <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  taxa_ref <- metagenomics$new(
    countData = taxa_hdf5$countData,
    metaData = taxa_hdf5$metaData,
    treeData = taxa_hdf5$treeData,
    featureData = taxa_hdf5$featureData
  )
  
  expect_snapshot(taxa_hdf5)
  expect_snapshot(taxa_ref)

    # Adding treeData after init
  taxa_ref <- metagenomics$new(
    countData = taxa_hdf5$countData,
    metaData = taxa_hdf5$metaData,
    featureData = taxa_hdf5$featureData,
  )

  expect_error(taxa_ref$treeData <- taxa_hdf5$countData)
  expect_no_error(taxa_ref$treeData <- taxa_hdf5$treeData)

  # Loading JSON format
  taxa_json <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_json.biom",
    metaData = data.table::data.table(SAMPLE_ID = c(
      "Sample1", "Sample2", "Sample3",
      "Sample4", "Sample5", "Sample6"))
  )
  expect_snapshot(taxa_json)
})
