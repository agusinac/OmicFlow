## Load example data
test <- omics$new(
  metaData = "input/proteomics/metadata.csv",
  countData = "input/proteomics/counts.csv"
)

test_that("`omics$scale()` -- Argument checks", {
  expect_snapshot(test$scale(method = "nonexistent"), error = TRUE)
  expect_snapshot(test$scale(method = 5), error = TRUE)
  expect_snapshot(test$scale(method = list()), error = TRUE)

  expect_snapshot(test$scale(transform = list()), error = TRUE)
  expect_snapshot(test$scale(transform = c(log2)), error = TRUE)

  expect_snapshot(test$scale(base = "1"), error = TRUE)
  expect_snapshot(test$scale(base = list()), error = TRUE)
  expect_snapshot(test$scale(base = c(1,2)), error = TRUE)

  expect_snapshot(test$scale(pseudocount = "1"), error = TRUE)
  expect_snapshot(test$scale(pseudocount = list()), error = TRUE)

  # This should not be allowed, but let's hope nobody will do such a thing
  # expect_snapshot(test$scale(pseudocount = c(1)), error = TRUE)
})

test_that("`omics$scale()` -- Behavioral checks", { 
  # Perform log transformation
  test$scale(method = "none", transform = log2)
  expect_snapshot(as.vector(test$countData[, 1]))
  
  # Perform sqrt transformation
  test$reset()
  test$scale(method = "none", transform = sqrt)
  expect_snapshot(as.vector(test$countData[, 1]))
  
  # Perform clr standardisation
  test$reset()
  mat <- as.matrix(test$countData)

  test$scale(method = "clr")
  expect_snapshot(as.vector(test$countData[, 1]))

  ## Comparing approaches
  ## vegan::decostand approach with `clr`
  clog <- log(mat)
  clog <- clog - rowMeans(clog, na.rm = TRUE)
  expect_equal(clog, as.matrix(test$countData))

  # Perform tss normalisation
  test$reset()
  test$scale(method = "tss")
  expect_snapshot(as.vector(test$countData[, 1]))

  # Perform hellinger transformation
  test$reset()
  test$scale(method = "hellinger")
  expect_snapshot(as.vector(test$countData[, 1]))

  # Perform binary transformation
  test$reset()
  test$scale(method = "binary")
  expect_snapshot(as.vector(test$countData[, 1]))

  # Perform clr normalisation with pseudocounts
  test$reset()
  test$scale(method = "clr", pseudocount = 1)
  expect_snapshot(as.vector(test$countData[, 1]))

  # Perform clr normalisation with different log base
  test$reset()
  test$scale(method = "clr", base = 2)
  expect_snapshot(as.vector(test$countData[, 1]))
})