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
  
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", method = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", method = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", aggregate_method = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", aggregate_method = "1"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", group_by = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", group_by = list()), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", group_by = "nothing"), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", group_by = c("1", "2")), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", feature_merge = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", feature_merge = "FALSE"), error = TRUE)

  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", paired = 1), error = TRUE)
  expect_snapshot(test$foldchange(condition.group = "treatment", condition_A = "tumor", condition_B = "healthy", paired = "FALSE"), error = TRUE)
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

    ## Adding `group_by`
    test$metaData[, "sex" := c("male", "male", "female", "male", "female", "male", "male", "female", "male", "female")]
    expect_no_error(res <- test$foldchange(condition.group = "CONTRAST_treatment", condition_A = "tumor", condition_B = "control", method = "log", group_by = "sex"))
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
})
