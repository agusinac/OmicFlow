test_that("Testing Proteomics loading", {
  counts_with_rownames <- proteomics$new(
    metaData = "input/proteomics/metadata.csv",
    countData = "input/proteomics/counts.csv",
    treeData = "input/proteomics/tree.newick"
  )

  counts_without_rownames <- proteomics$new(
    metaData = "input/proteomics/metadata.csv",
    countData = "input/proteomics/counts_without_rownames.csv"
  )
  
  prot_data <- proteomics$new(
    countData = counts_with_rownames$countData,
    metaData = counts_with_rownames$metaData,
    featureData = counts_with_rownames$featureData,
    treeData = counts_with_rownames$treeData
  )
  expect_snapshot(counts_without_rownames)
  expect_snapshot(counts_with_rownames)
  expect_snapshot(prot_data)

  # Adding treeData after init
  prot_1 <- proteomics$new(
    countData = counts_with_rownames$countData,
    metaData = counts_with_rownames$metaData,
    featureData = counts_with_rownames$featureData,
  )

  expect_error(prot_1$treeData <- counts_with_rownames$countData)
  expect_no_error(prot_1$treeData <- counts_with_rownames$treeData)
})
