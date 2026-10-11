## Load example data
test <- proteomics$new(
  metaData = "input/proteomics/metadata.csv",
  countData = "input/proteomics/counts.csv"
)

test_that("`hill_taxa()` -- Argument checks", {
  expect_snapshot(hill_taxa(x = data.frame()), error = TRUE)
  expect_snapshot(hill_taxa(x = c(2, 1, 1)), error = TRUE)
  expect_snapshot(hill_taxa(x = test$countData, normalize = "FALSE"), error = TRUE)
  expect_snapshot(hill_taxa(x = test$countData, base = "1"), error = TRUE)
  expect_snapshot(hill_taxa(x = test$countData, base = c(1, 2)), error = TRUE)
  expect_snapshot(hill_taxa(x = test$countData, q = "1"), error = TRUE)
  expect_snapshot(hill_taxa(x = test$countData, q = c(0, 1)), error = TRUE)

  ## Assuring values are positive
  test$countData[1,1] <- -1
  expect_snapshot(hill_taxa(x = test$countData), error = TRUE)
  test$reset()

})

test_that("`hill_taxa()` -- Behavioral checks", { 
  ## Testing different metrics
  # q = 0
  expect_no_error(res <- hill_taxa(x = test$countData, q = 0))
  expect_snapshot(res)

  # q = 1
  expect_no_error(res <- hill_taxa(x = test$countData, q = 1))
  expect_snapshot(res)

  # q = 2
  expect_no_error(res <- hill_taxa(x = test$countData, q = 2))
  expect_snapshot(res)

  ## Applying different base
  # base = 2
  expect_no_error(res <- hill_taxa(x = test$countData, q = 0, base = 2))
  expect_snapshot(res)

  # base = 10
  expect_no_error(res <- hill_taxa(x = test$countData, q = 0, base = 10))
  expect_snapshot(res)

  ## Applying with or no normalisation
  expect_no_error(res <- hill_taxa(x = test$countData, q = 0, normalize = FALSE))
  expect_snapshot(res)
})