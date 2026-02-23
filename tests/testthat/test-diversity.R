test_that("Testing `diversity` on sparse data", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  expect_no_error(diversity(taxa$countData, metric = "shannon"))
  expect_error(diversity(x = matrix_to_dtable(taxa$countData), metric = "shannon"))
  expect_error(diversity(taxa$countData, metric = 5))
  expect_error(diversity(taxa$countData, metric = "shan"))
})

test_that("Testing `diversity` on dense data", {
  prot <- proteomics$new(
      metaData = "input/proteomics/metadata.csv",
      countData = "input/proteomics/counts.csv"
  )
  
  expect_no_error(diversity(prot$countData, metric = "shannon"))
})