[![Codecov](https://codecov.io/gh/agusinac/OmicFlow/graph/badge.svg)](https://app.codecov.io/gh/agusinac/OmicFlow)
[![R-CMD-check](https://github.com/agusinac/OmicFlow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/agusinac/OmicFlow/actions/workflows/R-CMD-check.yaml)

OmicFlow
================

## Installation

The latest stable version can be installed from CRAN.

<!-- ``` r
install.packages('OmicFlow')
``` -->

The development version is available on GitHub.

``` r
install.packages('pak') # if not yet installed
pak::pak('agusinac/OmicFlow')
```

## Usage

### Simple loading and subsetting.
```R
library("OmicFlow")

taxa <- metagenomics$new(
    biomData = "tests/testthat/input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "tests/testthat/input/metagenomics/metadata.tsv",
    treeData = "tests/testthat/input/metagenomics/rooted_tree.newick"
)

taxa$feature_subset(Kingdom == "Bacteria")
taxa$normalize()

```

### Alpha- and Beta-diversity visualisations
```R
alpha_div <- taxa$alpha_diversity(
    col_name = "treatment",
    metric = "shannon"
)

alpha_div$diversity_plt$plot
```

```R
pcoa_plots <- taxa$ordination(
    metric = "unifrac",
    method = "pcoa",
    group_by = "treatment"
)

pcoa_plots
```

### Create an initial automated analysis in html format
```R
taxa$autoFlow(
    normalize = FALSE,
    weighted = TRUE,
    pvalue.threshold = 0.05
)
```

## Support

If you are having issues, please [create an issue](https://github.com/agusinac/OmicFlow/issues)
