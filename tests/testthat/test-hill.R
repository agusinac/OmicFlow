test_that("Testing Hill numbers on sparse data", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  res_0 <- hill_taxa(x = taxa$countData, q = 0)
  res_1 <- hill_taxa(x = taxa$countData, q = 1)
  res_2 <- hill_taxa(x = taxa$countData, q = 2)
  
  expect_snapshot(cat(res_0))
  expect_snapshot(cat(res_1))
  expect_snapshot(cat(res_2))
})

test_that("Testing Hill numbers on dense data", {
  prot <- proteomics$new(
      metaData = "input/proteomics/metadata.csv",
      countData = "input/proteomics/counts.csv"
  )
  
  res_0 <- hill_taxa(x = prot$countData, q = 0)
  res_1 <- hill_taxa(x = prot$countData, q = 1)
  res_2 <- hill_taxa(x = prot$countData, q = 2)
  
  expect_snapshot(cat(res_0))
  expect_snapshot(cat(res_1))
  expect_snapshot(cat(res_2))
})