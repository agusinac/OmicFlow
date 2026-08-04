## Load example data
test <- proteomics$new(
  metaData = "input/proteomics/metadata.csv",
  countData = "input/proteomics/counts.csv"
)

test_that("`diversity()` -- Argument checks", {
  expect_snapshot(diversity(x = data.frame()), error = TRUE)
  expect_snapshot(diversity(x = c(2, 1, 1)), error = TRUE)
  expect_snapshot(diversity(x = test$countData, normalize = "FALSE"), error = TRUE)
  expect_snapshot(diversity(x = test$countData, base = "1"), error = TRUE)
  expect_snapshot(diversity(x = test$countData, base = c(1, 2)), error = TRUE)
  expect_snapshot(diversity(x = test$countData, metric = 1), error = TRUE)
  expect_snapshot(diversity(x = test$countData, metric = c("shannon", "simpson")), error = TRUE)

  ## Assuring values are positive
  test$countData[1,1] <- -1
  expect_snapshot(diversity(x = test$countData), error = TRUE)
  test$reset()
})

test_that("`diversity()` -- Behavioral checks", { 
  ## Testing different metrics
  # Shannon
  expect_no_error(res <- diversity(x = test$countData, metric = "shannon"))
  expect_snapshot(res)

  # Simpson
  expect_no_error(res <- diversity(x = test$countData, metric = "simpson"))
  expect_snapshot(res)

  # Inverse Simpson
  expect_no_error(res <- diversity(x = test$countData, metric = "invsimpson"))
  expect_snapshot(res)

  ## Applying different base
  # base = 2
  expect_no_error(res <- diversity(x = test$countData, metric = "shannon", base = 2))
  expect_snapshot(res)

  # base = 10
  expect_no_error(res <- diversity(x = test$countData, metric = "shannon", base = 10))
  expect_snapshot(res)

  ## Applying with or no normalisation
  expect_no_error(res <- diversity(x = test$countData, metric = "shannon", normalize = FALSE))
  expect_snapshot(res)
})