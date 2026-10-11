## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`omics$autoFlow()` -- Argument checks", {
    expect_snapshot(taxa$autoFlow(feature_contrast = 1), error = TRUE)
    expect_snapshot(taxa$autoFlow(feature_contrast = c("2", "3")), error = TRUE)

    expect_snapshot(taxa$autoFlow(distmat = "path/nonexisting"), error = TRUE)
    expect_snapshot(taxa$autoFlow(distmat = matrix(0,0,0)), error = TRUE)

    expect_snapshot(taxa$autoFlow(report = "FALSE"), error = TRUE)
    expect_snapshot(taxa$autoFlow(filename = c("path/1", "path/2")), error = TRUE)

    colnames(taxa$metaData) <- c("SAMPLE_ID", "FW", "REV", "desc", "treat", "sex", "bmi", "weight", "time")
    expect_snapshot(taxa$autoFlow(), error = TRUE)
})

test_that("`omics$autoFlow()` -- Behavioral checks", { 
    taxa$reset()
    taxa$scale(method = "tss")
    expect_no_error(  
        suppressWarnings(
            taxa_autoflow <- taxa$autoFlow(
                feature_contrast = c("Phylum", "Family", "Genus"),
                pvalue.threshold = 1,
                report = FALSE
            )
    ))
    expect_true(length(taxa_autoflow$plots) == length(taxa_autoflow$data))
    expect_s3_class(taxa_autoflow$plots$alpha_div_plots[[1]], "ggplot")

    ## check feature contrasts
    expect_s3_class(taxa_autoflow$plots$composition_plots[[1]], "ggplot")
    expect_equal(ggplot2::ggplot_build(taxa_autoflow$plots$composition_plots[[1]])@plot@labels$fill, "Phylum")
    expect_s3_class(taxa_autoflow$plots$composition_plots[[2]], "ggplot")
    expect_equal(ggplot2::ggplot_build(taxa_autoflow$plots$composition_plots[[2]])@plot@labels$fill, "Family")
    expect_s3_class(taxa_autoflow$plots$composition_plots[[3]], "ggplot")
    expect_equal(ggplot2::ggplot_build(taxa_autoflow$plots$composition_plots[[3]])@plot@labels$fill, "Genus")

    expect_s3_class(taxa_autoflow$plots$pcoa_plots[[1]], "ggplot")
    expect_s3_class(taxa_autoflow$plots$Log2FC_plots[[1]], "ggplot")
    expect_s3_class(taxa_autoflow$plots$Log2FC_plots[[2]], "ggplot")
    expect_s3_class(taxa_autoflow$plots$Log2FC_plots[[3]], "ggplot")

    ## Run when skipping feature merge:
    expect_no_error(  
        suppressWarnings(
            taxa_autoflow <- taxa$autoFlow(
                feature_contrast = c("FEATURE_ID"),
                pvalue.threshold = 1,
                report = FALSE
            )
    ))
    expect_true(length(taxa_autoflow$plots) == length(taxa_autoflow$data))
    expect_s3_class(taxa_autoflow$plots$alpha_div_plots[[1]], "ggplot")

    ## check feature contrasts
    expect_s3_class(taxa_autoflow$plots$composition_plots[[1]], "ggplot")
    expect_equal(ggplot2::ggplot_build(taxa_autoflow$plots$composition_plots[[1]])@plot@labels$fill, "FEATURE_ID")
    expect_equal(length(taxa_autoflow$plots$composition_plots), 1)
    expect_equal(length(taxa_autoflow$plots$Log2FC_plots), 1)
    expect_s3_class(taxa_autoflow$plots$pcoa_plots[[1]], "ggplot")

    ## Run with custom distance matrix
    distmat <- taxa$distance(metric = "unifrac")

    expect_no_error(  
        suppressWarnings(
            taxa_autoflow <- taxa$autoFlow(
                feature_contrast = c("FEATURE_ID"),
                pvalue.threshold = 1,
                report = FALSE,
                distmat = as.matrix(distmat)
            )
    ))
    expect_true(length(taxa_autoflow$plots) == length(taxa_autoflow$data))
    expect_s3_class(taxa_autoflow$plots$alpha_div_plots[[1]], "ggplot")

    ## check feature contrasts
    expect_s3_class(taxa_autoflow$plots$composition_plots[[1]], "ggplot")
    expect_equal(ggplot2::ggplot_build(taxa_autoflow$plots$composition_plots[[1]])@plot@labels$fill, "FEATURE_ID")
    expect_equal(length(taxa_autoflow$plots$composition_plots), 1)
    expect_equal(length(taxa_autoflow$plots$Log2FC_plots), 1)

    expect_s3_class(taxa_autoflow$plots$pcoa_plots[[1]], "ggplot")
    expect_equal(taxa_autoflow$data$pcoa_data[[1]]$dist_mat[1,], as.matrix(distmat)[1,])
})