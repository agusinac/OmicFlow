## Load example data
test <- proteomics$new(
  metaData = "input/proteomics/metadata.csv",
  countData = "input/proteomics/counts.csv"
)
distmat <- test$distance(metric = "euclidean")

test_that("`pairwise_adonis()` -- Argument checks", {
  expect_snapshot(pairwise_adonis(x = c(1, 2, 3)), error = TRUE)
  expect_snapshot(pairwise_adonis(x = data.frame()), error = TRUE)
  expect_snapshot(pairwise_adonis(x = distmat, groups = list("4", "2")), error = TRUE)
  expect_snapshot(pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, metadata = matrix()), error = TRUE)
  expect_snapshot(pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, perm_design = matrix()), error = TRUE)
  expect_snapshot(pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, p.adjust.method = 1), error = TRUE)
  expect_snapshot(pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = c(1, 5)), error = TRUE)
  expect_snapshot(pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = "1"), error = TRUE)
  expect_snapshot(pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = 5.2), error = TRUE)
})

test_that("`pairwise_adonis()` -- Behavioral checks", {
  ## Testing without permutation design
  expect_no_error(res <- pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment))
  expect_equal(colnames(res), c("pairs", "Df", "SumsOfSqs", "F.Model", "R2", "p.value", "p.adj"))
  expect_equal(res$pairs, "tumor vs control")
  expect_equal(res$Df, 1)
  expect_equal(round(res$SumsOfSqs, 2), 1382.68)
  expect_equal(round(res$F.Model, 2), 1.68)
  expect_equal(round(res$R2, 2), 0.17)

  # ## Testing with permutation design
  # perm_design_func <- function(meta) {
  #   base::with(
  #     data = meta,
  #     expr = permute::how(
  #       nperm = 9,
  #       plots = permute::Plots(meta$SAMPLEPAIR_ID, type = "none"),
  #       within = permute::Within(type = "free")
  #     )
  #   )
  # }
  # expect_no_error(res <- pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, metadata = test$metaData, perm_design = perm_design_func))
  # expect_equal(colnames(res), c("pairs", "Df", "SumsOfSqs", "F.Model", "R2", "p.value", "p.adj"))
  # expect_equal(res$pairs, "tumor vs control")
  # expect_equal(res$Df, 1)
  # expect_equal(round(res$SumsOfSqs, 2), 1382.68)
  # expect_equal(round(res$F.Model, 2), 1.68)
  # expect_equal(round(res$R2, 2), 0.17)
  # expect_equal(res$p.value, 1)
})

test_that("`pairwise_anosim()` -- Argument checks", {
  expect_snapshot(pairwise_anosim(x = c(1, 2, 3)), error = TRUE)
  expect_snapshot(pairwise_anosim(x = data.frame()), error = TRUE)
  expect_snapshot(pairwise_anosim(x = distmat, groups = list("4", "2")), error = TRUE)
  expect_snapshot(pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, metadata = matrix()), error = TRUE)
  expect_snapshot(pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, perm_design = matrix()), error = TRUE)
  expect_snapshot(pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, p.adjust.method = 1), error = TRUE)
  expect_snapshot(pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = c(1, 5)), error = TRUE)
  expect_snapshot(pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = "1"), error = TRUE)
  expect_snapshot(pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = 5.2), error = TRUE)
})

test_that("`pairwise_anosim()` -- Behavioral checks", { 
  ## Testing without permutation design
  expect_no_error(res <- pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment))
  expect_equal(colnames(res), c("pairs", "anosimR", "p.value", "p.adj"))
  expect_equal(res$pairs, "tumor vs control")
  expect_equal(round(res$anosimR, 2), 0.38)

  # ## Testing with permutation design
  # perm_design_func <- function(meta) {
  #   base::with(
  #     data = meta,
  #     expr = permute::how(
  #       nperm = 9,
  #       plots = permute::Plots(meta$SAMPLEPAIR_ID, type = "none"),
  #       within = permute::Within(type = "free")
  #     )
  #   )
  # }
  # expect_no_error(res <- pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, metadata = test$metaData, perm_design = perm_design_func))
  # expect_equal(colnames(res), c("pairs", "anosimR", "p.value", "p.adj"))
  # expect_equal(res$pairs, "tumor vs control")
  # expect_equal(round(res$anosimR, 2), 0.38)
  # expect_equal(res$p.value, 1)
})


dt <- test$alpha_diversity(col_name = "CONTRAST_treatment")$data

test_that("`pairwise_wilcox_test()` -- Argument checks", {
  expect_snapshot(pairwise_wilcox_test(data = list()), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = matrix()), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = 1), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = c("1", "2")), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "nonexisting"), error = TRUE)

  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "V1", g_col = 1), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "V1", g_col = c("1", "2")), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "nonexisting"), error = TRUE)

  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "CONTRAST_treatment", paired = 1), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "CONTRAST_treatment", paired = "FALSE"), error = TRUE)

  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "CONTRAST_treatment", p.adjust.method = 1), error = TRUE)
  expect_snapshot(pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "CONTRAST_treatment", p.adjust.method = "nonexisting"), error = TRUE)
})

test_that("`pairwise_wilcox_test()` -- Behavioral checks", {
  ## Testing with default settings
  expect_no_error(res <- pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "CONTRAST_treatment"))
  expect_equal(res$group1, "tumor")
  expect_equal(res$group2, "control")
  expect_equal(colnames(res)[5], "obs.tot")
  expect_equal(res$statistic, 12)
  expect_equal(res$pvalue, 1)
  expect_equal(res$xmin, 1)
  expect_equal(res$xmax, 2)

  ## Testing with `paired = TRUE`
  expect_no_error(res <- pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "CONTRAST_treatment", paired = TRUE))
  expect_equal(res$group1, "tumor")
  expect_equal(res$group2, "control")
  expect_equal(colnames(res)[5], "obs.paired")
  expect_equal(res$statistic, 5)
  expect_equal(res$pvalue, 0.625)
  expect_equal(res$xmin, 1)
  expect_equal(res$xmax, 2)

  ## Testing with multiple groups
  dt[[ "CONTRAST_treatment" ]][4:5] <- c("unknown")
  dt[[ "CONTRAST_treatment" ]][9:10] <- c("unknown")
  expect_no_error(res <- pairwise_wilcox_test(data = dt, x_col = "V1", g_col = "CONTRAST_treatment"))
  expect_equal(res$group1, c("tumor", "tumor", "unknown"))
  expect_equal(res$group2, c("unknown", "control", "control"))
  expect_equal(colnames(res)[5], "obs.tot")
  expect_equal(res$statistic, c(8, 4, 2))
  expect_equal(round(res$pvalue, 2), c(0.63, 1.00, 0.23))
  expect_equal(res$xmin, c(1, 1, 2))
  expect_equal(res$xmax, c(3, 2, 2))
})