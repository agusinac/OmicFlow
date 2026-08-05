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
  ## Testing default setting with method 'pcoa'
  expect_no_error(res <- test$ordination(group_by = "CONTRAST_treatment"))

  expect_equal(res$anova_data$pairs, "tumor vs control")
  expect_equal(round(res$anova_data$F.Model, 2), 1.75)
  expect_equal(dim(as.matrix(res$dist)), c(nrow(test$metaData), nrow(test$metaData)))
  expect_equal(any(grepl("PC", colnames(res$pcs))), TRUE)
  expect_equal(column_exists("groups", res$pcs), TRUE)
  expect_equal(column_exists("samples", res$pcs), TRUE)

  expect_s3_class(res$scores_plot, "ggplot")
  expect_s3_class(res$scree_plot, "ggplot")
  expect_s3_class(res$anova_plot, "ggplot")

  ## Testing setting with method 'nmds'
  expect_no_error(res <- suppressWarnings(test$ordination(group_by = "CONTRAST_treatment", method = "nmds")))

  expect_equal(res$anova_data$pairs, "tumor vs control")
  expect_equal(round(res$anova_data$anosimR, 2), 0.32)
  expect_equal(dim(as.matrix(res$dist)), c(nrow(test$metaData), nrow(test$metaData)))
  expect_equal(any(grepl("MDS", colnames(res$pcs))), TRUE)
  expect_equal(column_exists("stress", res$pcs), TRUE)
  expect_equal(column_exists("groups", res$pcs), TRUE)
  expect_equal(column_exists("samples", res$pcs), TRUE)

  expect_s3_class(res$scores_plot, "ggplot")
  expect_s3_class(res$anova_plot, "ggplot")
  
  ## Check if distmat can be supplied as `Matrix` or `dist` class
  distmat <- test$distance(metric = "canberra")

  expect_no_error(res <- test$ordination(group_by = "CONTRAST_treatment", distmat = distmat))

  expect_equal(res$anova_data$pairs, "tumor vs control")
  expect_equal(round(res$anova_data$F.Model, 2), 1.74)
  expect_equal(dim(as.matrix(res$dist)), c(nrow(test$metaData), nrow(test$metaData)))
  expect_equal(any(grepl("PC", colnames(res$pcs))), TRUE)
  expect_equal(column_exists("groups", res$pcs), TRUE)
  expect_equal(column_exists("samples", res$pcs), TRUE)

  expect_s3_class(res$scores_plot, "ggplot")
  expect_s3_class(res$scree_plot, "ggplot")
  expect_s3_class(res$anova_plot, "ggplot")
})