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


## Load example data
n_row <- 1000
n_col <- 100
density <- 0.2
num_entries <- n_row * n_col
num_nonzero <- round(num_entries * density)

set.seed(123)
positions <- sample(num_entries, num_nonzero, replace=FALSE)
row_idx <- ((positions - 1) %% n_row) + 1
col_idx <- ((positions - 1) %/% n_row) + 1

values <- runif(num_nonzero, min = 0, max = 1)
sparse_mat <- Matrix::sparseMatrix(
  i = row_idx,
  j = col_idx,
  x = values,
  dims = c(n_row, n_col)
)

div <- OmicFlow::diversity(
x = sparse_mat,
metric = "shannon"
)

dt <- data.table::data.table(
"shannon" = div,
"treatment" = c(rep("healthy", n_col / 2), rep("tumor", n_col / 2)),
"sex" = c(rep("male", n_col / 4), rep("female", n_col / 4))
)

test_that("`pairwise_test()` -- Argument checks", {
  expect_snapshot(pairwise_test(data = list()), error = TRUE)
  expect_snapshot(pairwise_test(data = matrix()), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = 1), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = c("1", "2")), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = "nonexisting"), error = TRUE)

  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = 1), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = c("1", "2")), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = "nonexisting"), error = TRUE)

  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", paired = 1), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", paired = "FALSE"), error = TRUE)

  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", p.adjust.method = 1), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", p.adjust.method = "nonexisting"), error = TRUE)

  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = 1), error = TRUE)
  expect_snapshot(pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "nonexisting"), error = TRUE)
})

test_that("`pairwise_wilcox_test()` -- Behavioral checks", {
  ## Testing with default settings
  expect_no_error(res <- pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "wilcox"))
  expect_equal(res$group1, "healthy")
  expect_equal(res$group2, "tumor")
  expect_equal(colnames(res)[5], "obs.tot")
  expect_equal(res$statistic, 1453)
  expect_equal(round(res$pvalue, 2), 0.16)
  expect_equal(res$xmin, 1)
  expect_equal(res$xmax, 2)

  ## Testing with `paired = TRUE`
  expect_no_error(res <- pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "wilcox", paired = TRUE))
  expect_equal(res$group1, "healthy")
  expect_equal(res$group2, "tumor")
  expect_equal(colnames(res)[5], "obs.paired")
  expect_equal(res$statistic, 749)
  expect_equal(round(res$pvalue, 2), 0.28)
  expect_equal(res$xmin, 1)
  expect_equal(res$xmax, 2)

  
  ## Testing t.test paired
  expect_no_error(res <- pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "t.test", paired = TRUE))
  expect_equal(res$group1, "healthy")
  expect_equal(res$group2, "tumor")
  expect_equal(colnames(res)[5], "obs.paired")
  expect_equal(round(res$statistic, 2), 1.13)
  expect_equal(round(res$mean.diff, 3), 0.015)
  expect_equal(round(res$var.diff, 3), 0.009)
  expect_equal(round(res$stderr, 3), 0.013)
  expect_equal(res$xmin, 1)
  expect_equal(res$xmax, 2)

  ## Testing t.test unpaired
  expect_no_error(res <- pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "t.test", paired = FALSE))
  expect_equal(res$group1, "healthy")
  expect_equal(res$group2, "tumor")
  expect_equal(colnames(res)[5], "obs.tot")
  expect_equal(round(res$statistic, 2), 1.24)
  expect_equal(round(res$mean.diff, 3), 0.015)
  expect_equal(round(res$stderr, 3), 0.012)
  expect_equal(res$xmin, 1)
  expect_equal(res$xmax, 2)

  ## Testing t.equalvar unpaired
  expect_no_error(res <- pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "t.equalvar", paired = FALSE))
  expect_equal(res$group1, "healthy")
  expect_equal(res$group2, "tumor")
  expect_equal(colnames(res)[5], "obs.tot")
  expect_equal(round(res$statistic, 2), 1.24)
  expect_equal(round(res$mean.diff, 3), 0.015)
  expect_equal(round(res$var.pooled, 3), 0.004)
  expect_equal(round(res$stderr, 3), 0.012)
  expect_equal(res$xmin, 1)
  expect_equal(res$xmax, 2)

  ## Testing with multiple groups
  dt[[ "treatment" ]][1:10] <- c("unknown")
  dt[[ "treatment" ]][90:100] <- c("unknown")
  expect_no_error(res <- pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "wilcox"))
  expect_equal(res$group1, c("unknown", "unknown", "healthy"))
  expect_equal(res$group2, c("healthy", "tumor", "tumor"))
  expect_equal(colnames(res)[5], "obs.tot")
  expect_equal(res$statistic, c(422, 465, 908))
  expect_equal(round(res$pvalue, 2), c(0.98, 0.40, 0.21))
  expect_equal(res$xmin, c(2, 2, 1))
  expect_equal(res$xmax, c(2, 3, 3))
})