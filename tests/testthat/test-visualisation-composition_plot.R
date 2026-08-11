## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)
taxa$scale(method = "tss")

comp <- taxa$composition(
  feature_rank = "Genus",
  feature_filter = c("uncultured"),
  col_name = "CONTRAST_sex",
  feature_top = 10
)

test_that("`composition_plot()` -- Argument checks", {
  expect_snapshot(composition_plot(data = list()), error = TRUE)
  expect_snapshot(composition_plot(data = matrix()), error = TRUE)
  expect_snapshot(composition_plot(data = comp$data, palette = list()), error = TRUE)
  expect_snapshot(composition_plot(data = comp$data, palette = c("foo1", "foo2")), error = TRUE)
  expect_snapshot(composition_plot(data = comp$data, palette = comp$palette, feature_rank = 1), error = TRUE)
  expect_snapshot(composition_plot(data = comp$data, palette = comp$palette, feature_rank = c("1", "2")), error = TRUE)

  expect_snapshot(composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus", title_name = 1), error = TRUE)
  expect_snapshot(composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus", title_name = c("1", "2")), error = TRUE)

  expect_snapshot(composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus", group_by = 1), error = TRUE)
  expect_snapshot(composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus", group_by = c("1", "2")), error = TRUE)
})


test_that("`composition_plot()` -- Behavioral checks", {
  ## Testing default
  expect_no_error(p <- composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus"))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@data[1, ], comp$data[1, ])
  expect_equal(res@plot@labels$y, "Rel. Abun.")

  ## Changing `title_name`
  title_name <- "New title"
  expect_no_error(p <- composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus", title_name = title_name))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@data[1, ], comp$data[1, ])
  expect_equal(res@plot@labels$y, "Rel. Abun.")
  expect_equal(res@plot@labels$title, title_name)

  ## Changing `group_by`
  group_by <- "CONTRAST_sex"
  expect_no_error(p <- composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus", group_by = group_by))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@data[1, ], comp$data[1, ])
  expect_equal(res@plot@labels$y, "Rel. Abun.")
})