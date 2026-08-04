# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library("testthat")
library("OmicFlow")
set.seed(100)

## Load example data
metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_sparse_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
counts_sparse_with_rownames_file <- system.file("extdata", "counts_with_rownames.tsv", package = "OmicFlow")
counts_dense_file <- system.file("extdata", "counts_dense.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

## Start tests per section
test_check("OmicFlow", filter = "omics")
test_check("OmicFlow", filter = "metagenomics")
test_check("OmicFlow", filter = "proteomics")
test_check("OmicFlow", filter = "metrics")
test_check("OmicFlow", filter = "statistics")
# test_check("OmicFlow", filter = "visualisation")
test_check("OmicFlow", filter = "util")
