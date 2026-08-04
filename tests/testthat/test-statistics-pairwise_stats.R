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

  ## Testing with permutation design
  perm_design_func <- function(meta) {
    base::with(
      data = meta,
      expr = permute::how(
        nperm = 9,
        plots = permute::Plots(meta$SAMPLEPAIR_ID, type = "none"),
        within = permute::Within(type = "free")
      )
    )
  }
  expect_no_error(res <- pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, metadata = test$metaData, perm_design = perm_design_func))
  expect_equal(colnames(res), c("pairs", "Df", "SumsOfSqs", "F.Model", "R2", "p.value", "p.adj"))
  expect_equal(res$pairs, "tumor vs control")
  expect_equal(res$Df, 1)
  expect_equal(round(res$SumsOfSqs, 2), 1382.68)
  expect_equal(round(res$F.Model, 2), 1.68)
  expect_equal(round(res$R2, 2), 0.17)
  expect_equal(res$p.value, 1)
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

  ## Testing with permutation design
  perm_design_func <- function(meta) {
    base::with(
      data = meta,
      expr = permute::how(
        nperm = 9,
        plots = permute::Plots(meta$SAMPLEPAIR_ID, type = "none"),
        within = permute::Within(type = "free")
      )
    )
  }
  expect_no_error(res <- pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, metadata = test$metaData, perm_design = perm_design_func))
  expect_equal(colnames(res), c("pairs", "anosimR", "p.value", "p.adj"))
  expect_equal(res$pairs, "tumor vs control")
  expect_equal(round(res$anosimR, 2), 0.38)
  expect_equal(res$p.value, 1)
})