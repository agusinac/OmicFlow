## Load example data
test <- proteomics$new(
  metaData = "input/proteomics/metadata.csv",
  countData = "input/proteomics/counts.csv"
)

test_that("`omics$distance()` -- Argument checks", {
  expect_snapshot(test$distance(), error = TRUE)
  expect_snapshot(test$distance(metric = "nothing"), error = TRUE)
  expect_snapshot(test$distance(metric = c("n1", "n2")), error = TRUE)
  expect_snapshot(test$distance(metric = 1), error = TRUE)
  expect_snapshot(test$distance(metric = "unifrac"), error = TRUE)
  expect_snapshot(test$distance(metric = "bray", threads = "1"), error = TRUE)
  expect_snapshot(test$distance(metric = "bray", threads = 50.2), error = TRUE)
})

test_that("`omics$ordination()` -- Argument checks", {
  expect_snapshot(test$ordination(), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", metric = "nothing"), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", metric = c("n1", "n2")), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", metric = 1), error = TRUE)

  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", method = "nothing"), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", method = c("n1", "n2")), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", method = 1), error = TRUE)

  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", perm = "999"), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", perm = 50.2), error = TRUE)

  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", threads = "999"), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", threads = 50.2), error = TRUE)

  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", perm_design = list()), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", perm_design = function(x) print(x)), error = TRUE)

  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", distmat = list()), error = TRUE)
  expect_snapshot(test$ordination(group_by = "CONTRAST_treatment", distmat = dist(c(2,1,2))), error = TRUE)
})

test_that("`omics$ordination()` -- Behavioral checks", { 
  expect_no_error(res <- test$ordination(group_by = "CONTRAST_treatment"))

  expect_snapshot(res$anova_data)
  expect_snapshot(res$dist)
  expect_snapshot(res$pcs)

  expect_s3_class(res$scores_plot, "ggplot")
  expect_s3_class(res$scree_plot, "ggplot")
  expect_s3_class(res$anova_plot, "ggplot")
  
  ## Check if distmat can be supplied as `Matrix` or `dist` class
  distmat <- test$distance(metric = "canberra")

  expect_no_error(res <- test$ordination(group_by = "CONTRAST_treatment", distmat = distmat))

  expect_snapshot(res$anova_data)
  expect_snapshot(res$dist)
  expect_snapshot(res$pcs)

  expect_s3_class(res$scores_plot, "ggplot")
  expect_s3_class(res$scree_plot, "ggplot")
  expect_s3_class(res$anova_plot, "ggplot")
})