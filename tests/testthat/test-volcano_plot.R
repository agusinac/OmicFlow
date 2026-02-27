test_that("Testing `volcano_plot`", {
# Create mock data frame
mock_volcano_data <- data.table(

  # Feature names (feature_rank)
  Feature = paste0("Gene", 1:20),

  # Log2 fold changes (X)        
  log2FC = c(1.2, -1.5, 0.3, -0.7, 2.3,
             -2.0, 0.1, 0.5, -1.0, 1.8,
             -0.4, 0.7, -1.4, 1.5, 0.9,
             -2.1, 0.2, 1.0, -0.3, -1.8),
  
  # P-values (Y)
  pvalue = c(0.001, 0.02, 0.3, 0.04, 0.0005,
             0.01, 0.7, 0.5, 0.02, 0.0008,
             0.15, 0.06, 0.01, 0.005, 0.3,
             0.02, 0.8, 0.04, 0.12, 0.03),
  
  # Mean (relative) abundance for point sizing
  rel_abun = runif(20, 0.01, 0.1)
)

  expect_error(volcano_plot(
    data = as.matrix(mock_volcano_data),
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "val2",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pval",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "abun",
    feature_rank = "Feature",
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feat",
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
    pvalue.threshold = "low"
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
    logfold.threshold = "low"
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
    abundance.threshold = "low"
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
    label_A = 1
  ))

  expect_error(volcano_plot(
    data = mock_volcano_data,
    logfold_col = "log2FC",
    pvalue_col = "pvalue",
    abundance_col = "rel_abun",
    feature_rank = "Feature",
    label_B = 1
  ))
})