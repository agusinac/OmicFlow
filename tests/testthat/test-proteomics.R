test_that("Testing Proteomics loading", {
  prot_file <- proteomics$new(
    metaData = "input/proteomics/metadata.csv",
    countData = "input/proteomics/counts.csv",
    treeData = "input/proteomics/tree.newick"
  )
  
  prot_data <- proteomics$new(
    countData = prot_file$countData,
    metaData = prot_file$metaData,
    featureData = prot_file$featureData,
    treeData = prot_file$treeData
  )
  
  expect_snapshot(prot_file)
  expect_snapshot(prot_data)

  # Adding treeData after init
  prot_1 <- proteomics$new(
    countData = prot_file$countData,
    metaData = prot_file$metaData,
    featureData = prot_file$featureData,
  )

  expect_error(prot_1$treeData <- prot_file$countData)
  expect_no_error(prot_1$treeData <- prot_file$treeData)
})
