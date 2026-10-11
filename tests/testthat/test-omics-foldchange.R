## Load example data
test <- proteomics$new(
  metaData = "input/proteomics/metadata.csv",
  countData = "input/proteomics/counts.csv"
)

test_that("`omics$foldchange()` -- Argument checks", {
  expect_error(test$foldchange(), 'argument "condition.group" is missing, with no default')
  expect_snapshot(test$foldchange(condition.group = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = list()), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "nothing"), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = c("1", "2")), error = TRUE)
  
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", aggregate_method = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", aggregate_method = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", split_by = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", split_by = list()), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", split_by = "nothing"), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", split_by = c("1", "2")), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", feature_merge = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", feature_merge = "FALSE"), error = TRUE)

  ## TODO: `testthat` can't capture a `cli` warning message..
  # expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", paired = 1), error = TRUE)
  # expect_snapshot(test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", paired = "FALSE"), error = TRUE)
})

test_that("`omics$foldchange()` -- Behavioral checks", { 
  ## Testing with default settings
  ### Using `method = "identity"`
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control"))
  expect_s3_class(res$volcano_plot[[1]], "ggplot")
  ## returned values should be identical
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor"]$SAMPLE_ID]), res$all_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control"]$SAMPLE_ID]), res$all_control)
  expect_equal(colnames(res$data), c("FEATURE_ID", "median_abun", "fold-change_tumor_vs_control_in_all", "pvalue_wilcox_tumor_vs_control_in_all", "homogeneity_test_statistic_tumor_vs_control_in_all", "homogeneity_test_pvalue_tumor_vs_control_in_all"))
  expect_snapshot(res$data)

  ### Using `method = "identity"` with proportions
  test$scale(method = "tss")
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control"))
  expect_s3_class(res$volcano_plot[[1]], "ggplot")
  ## returned values should be identical
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor"]$SAMPLE_ID]), res$all_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control"]$SAMPLE_ID]), res$all_control)
  expect_equal(colnames(res$data), c("FEATURE_ID", "median_abun", "fold-change_tumor_vs_control_in_all", "pvalue_wilcox_tumor_vs_control_in_all", "homogeneity_test_statistic_tumor_vs_control_in_all", "homogeneity_test_pvalue_tumor_vs_control_in_all"))
  expect_snapshot(res$data)

  ### Using `method = "log"`
  test$reset()
  test$scale(method = "clr")
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "log"))
  expect_s3_class(res$volcano_plot[[1]], "ggplot")
  ## returned values should be identical
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor"]$SAMPLE_ID]), res$all_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control"]$SAMPLE_ID]), res$all_control)
  expect_equal(colnames(res$data), c("FEATURE_ID", "median_abun", "fold-change_tumor_vs_control_in_all", "pvalue_wilcox_tumor_vs_control_in_all", "homogeneity_test_statistic_tumor_vs_control_in_all", "homogeneity_test_pvalue_tumor_vs_control_in_all"))
  expect_snapshot(res$data)

  ## Adding `split_by`
  test$metaData[, "sex" := c("male", "male", "female", "male", "female", "male", "male", "female", "male", "female")]
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "log", split_by = "sex"))
  expect_s3_class(res$volcano_plot[[1]], "ggplot")
  expect_s3_class(res$volcano_plot[[2]], "ggplot")
  ## returned values should be identical
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor" & sex == "male"]$SAMPLE_ID]), res$male_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor" & sex == "female"]$SAMPLE_ID]), res$female_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control" & sex == "male"]$SAMPLE_ID]), res$male_control)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control" & sex == "female"]$SAMPLE_ID]), res$female_control)
  expect_snapshot(res$data)

  ## Setting `paired = TRUE`
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "log", paired = TRUE))
  expect_s3_class(res$volcano_plot[[1]], "ggplot")
  ## returned values should be identical
  expect_equal(colnames(res$data), c("FEATURE_ID", "median_abun", "fold-change_tumor_vs_control_in_all", "pvalue_wilcox-paired_tumor_vs_control_in_all", "homogeneity_test_statistic_tumor_vs_control_in_all", "homogeneity_test_pvalue_tumor_vs_control_in_all"))
  expect_snapshot(res$data)

  ## Setting `aggregate_method` to `mean`
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "log", aggregate_method = "mean"))
  expect_s3_class(res$volcano_plot[[1]], "ggplot")
  ## returned values should be identical
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor"]$SAMPLE_ID]), res$all_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control"]$SAMPLE_ID]), res$all_control)
  expect_equal(colnames(res$data), c("FEATURE_ID", "median_abun", "fold-change_tumor_vs_control_in_all", "pvalue_wilcox_tumor_vs_control_in_all", "homogeneity_test_statistic_tumor_vs_control_in_all", "homogeneity_test_pvalue_tumor_vs_control_in_all"))
  expect_snapshot(res$data)

  ## Setting `aggregate_method` to `geomean`, similar to `DESeq2`
  test$reset()
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "log", aggregate_method = "geomean"))
  expect_s3_class(res$volcano_plot[[1]], "ggplot")
  ## returned values should be identical
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor"]$SAMPLE_ID]), res$all_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control"]$SAMPLE_ID]), res$all_control)
  expect_equal(colnames(res$data), c("FEATURE_ID", "median_abun", "fold-change_tumor_vs_control_in_all", "pvalue_wilcox_tumor_vs_control_in_all", "homogeneity_test_statistic_tumor_vs_control_in_all", "homogeneity_test_pvalue_tumor_vs_control_in_all"))
  expect_snapshot(res$data)

  ## Setting `aggregate_method` to `none` (returns sample ids)
  expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "identity", aggregate_method = "none"))
  expect_equal(names(res), c("all_tumor", "all_control", "data"))
  ## returned values should be identical
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "tumor"]$SAMPLE_ID]), res$all_tumor)
  expect_equal(as.matrix(test$countData[, test$metaData[CONTRAST_treatment == "control"]$SAMPLE_ID]), res$all_control)
  sample_ids_A <- test$metaData[CONTRAST_treatment == "tumor"]$SAMPLE_ID

  ## Since `aggregate_method == "none"` homogeneity test is not applied since no aggregation method is used..
  ## This might be changed in the future if a specific issue is raised.
  expect_equal(colnames(res$data), c("FEATURE_ID", "median_abun", sample_ids_A, "pvalue_wilcox_tumor_vs_control_in_all"))
  expect_snapshot(res$data)
})