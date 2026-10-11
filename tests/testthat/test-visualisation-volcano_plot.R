## Load example data
mock_data <- data.table::data.table(
  Feature = paste0("Gene", 1:20),
  log2FC = c(1.2, -1.5, 0.3, -0.7, 2.3,
              -2.0, 0.1, 0.5, -1.0, 1.8,
              -0.4, 0.7, -1.4, 1.5, 0.9,
              -2.1, 0.2, 1.0, -0.3, -1.8),
  pvalue = c(0.001, 0.02, 0.3, 0.04, 0.0005,
              0.01, 0.7, 0.5, 0.02, 0.0008,
              0.15, 0.06, 0.01, 0.005, 0.3,
              0.02, 0.8, 0.04, 0.12, 0.03),
  rel_abun = runif(20, 0.01, 0.1)
)

test_that("`volcano_plot()` -- Argument checks", {
  expect_snapshot(volcano_plot(data = data.frame()), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = 1), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = c("1", "2")), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "nothing"), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = 1), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = c("1", "2")), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "nothing"), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = 1), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = c("1", "2")), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "nothing"), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = 1), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = c("1", "2")), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "nothing"), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", pvalue.threshold = "1"), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", pvalue.threshold = c(1, 2)), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", logfold.threshold = "1"), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", logfold.threshold = c(1, 2)), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", abundance.threshold = "1"), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", abundance.threshold = c(1, 2)), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", label_A = 1), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", label_A = c("A1", "A2")), error = TRUE)

  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", label_B = 1), error = TRUE)
  expect_snapshot(volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", label_B = c("A1", "A2")), error = TRUE)
})

test_that("`volcano_plot()` -- Behavioral checks", {
  ## Testing default
  expect_no_error(p <- volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature"))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(colnames(res@plot@data), c("Feature", "log2FC", "pvalue", "rel_abun", "diffexpressed", "diffexpressed_labels"))
  expect_equal(mock_data[log2FC > 0.6 & pvalue < 0.05]$Feature, res@plot@data[diffexpressed == "Upregulated"]$Feature)
  expect_equal(mock_data[log2FC < -0.6 & pvalue < 0.05]$Feature, res@plot@data[diffexpressed == "Downregulated"]$Feature)
  expect_equal(mock_data[pvalue > 0.05]$Feature, res@plot@data[diffexpressed == "non-significant"]$Feature)

  expect_equal(res@plot@labels$x, paste0("Fold Change ( A / B )"))
  expect_equal(res@plot@labels$y, paste0("-log10( pvalue )"))
  expect_equal(res@plot@labels$colour, paste0("log2FC"))
  expect_equal(res@plot@labels$label, paste0("diffexpressed_labels"))

  ## Changing `pvalue.threshold`
  pvalue_threshold <- 0.01
  expect_no_error(p <- volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", pvalue.threshold = pvalue_threshold))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(colnames(res@plot@data), c("Feature", "log2FC", "pvalue", "rel_abun", "diffexpressed", "diffexpressed_labels"))
  expect_equal(mock_data[log2FC > 0.6 & pvalue < pvalue_threshold]$Feature, res@plot@data[diffexpressed == "Upregulated"]$Feature)
  expect_equal(mock_data[log2FC < -0.6 & pvalue < pvalue_threshold]$Feature, res@plot@data[diffexpressed == "Downregulated"]$Feature)
  expect_equal(mock_data[pvalue >= pvalue_threshold]$Feature, res@plot@data[diffexpressed == "non-significant"]$Feature)

  expect_equal(res@plot@labels$x, paste0("Fold Change ( A / B )"))
  expect_equal(res@plot@labels$y, paste0("-log10( pvalue )"))
  expect_equal(res@plot@labels$colour, paste0("log2FC"))
  expect_equal(res@plot@labels$label, paste0("diffexpressed_labels"))

  ## Changing `logfold.threshold`
  logfold_threshold <- 0.1
  expect_no_error(p <- volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", logfold.threshold = logfold_threshold))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(colnames(res@plot@data), c("Feature", "log2FC", "pvalue", "rel_abun", "diffexpressed", "diffexpressed_labels"))
  expect_equal(mock_data[log2FC > logfold_threshold & pvalue < 0.05]$Feature, res@plot@data[diffexpressed == "Upregulated"]$Feature)
  expect_equal(mock_data[log2FC < -logfold_threshold & pvalue < 0.05]$Feature, res@plot@data[diffexpressed == "Downregulated"]$Feature)
  expect_equal(mock_data[pvalue > 0.05]$Feature, res@plot@data[diffexpressed == "non-significant"]$Feature)

  expect_equal(res@plot@labels$x, paste0("Fold Change ( A / B )"))
  expect_equal(res@plot@labels$y, paste0("-log10( pvalue )"))
  expect_equal(res@plot@labels$colour, paste0("log2FC"))
  expect_equal(res@plot@labels$label, paste0("diffexpressed_labels"))

  ## Changing `label_A` and `label_B`
  label_A <- "T1"
  label_B <- "T2"
  expect_no_error(p <- volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue", abundance_col = "rel_abun", feature_rank = "Feature", label_A = label_A, label_B = label_B))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(colnames(res@plot@data), c("Feature", "log2FC", "pvalue", "rel_abun", "diffexpressed", "diffexpressed_labels"))
  expect_equal(mock_data[log2FC > 0.6 & pvalue < 0.05]$Feature, res@plot@data[diffexpressed == "Upregulated"]$Feature)
  expect_equal(mock_data[log2FC < -0.6 & pvalue < 0.05]$Feature, res@plot@data[diffexpressed == "Downregulated"]$Feature)
  expect_equal(mock_data[pvalue > 0.05]$Feature, res@plot@data[diffexpressed == "non-significant"]$Feature)

  expect_equal(res@plot@labels$x, paste0("Fold Change ( ", label_A," / ", label_B, " )"))
  expect_equal(res@plot@labels$y, paste0("-log10( pvalue )"))
  expect_equal(res@plot@labels$colour, paste0("log2FC"))
  expect_equal(res@plot@labels$label, paste0("diffexpressed_labels"))
})