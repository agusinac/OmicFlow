# Ordination plot

Creates an ordination plot from pre-computed principal components. This
function is built into the class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) with
method `ordination()` and inherited by other omics classes, such as;
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md).

## Usage

``` r
ordination_plot(
  data,
  col_name = lifecycle::deprecated(),
  groups = lifecycle::deprecated(),
  pair,
  dist_explained = NULL,
  dist_metric = NULL
)
```

## Arguments

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) or
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) with
  Principal Components as columns and rows as loading scores.

- col_name:

  **\[deprecated\]** This argument has been renamed to `groups` for more
  clarity.

- groups:

  A categorical variable in `data` to color the groups (e.g.
  `"treatment"`).

- pair:

  A vector of character variables indicating what dimension names (e.g.
  `c("PC1", "PC2")` or `c("NMDS1", "NMDS2")`) in `data`.

- dist_explained:

  A vector of numeric values of the percentage dissimilarity explained
  for the dimension pairs (default: `NULL`).

- dist_metric:

  A character variable indicating what metric is used (e.g.
  `"bray-curtis"`) (default: `NULL`).

## Value

A
[ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object to be further modified

## Examples

``` r
# Mock principal component scores
set.seed(123)
mock_data <- data.frame(
  SampleID = paste0("Sample", 1:10),
  PC1 = rnorm(10, mean = 0, sd = 1),
  PC2 = rnorm(10, mean = 0, sd = 1),
  groups = rep(c("Group1", "Group2"), each = 5)
)

# Basic usage
ordination_plot(
  data = mock_data,
  groups = "groups",
  pair = c("PC1", "PC2")
)


# Adding variance/dissimilarity explained.
ordination_plot(
  data = mock_data,
  groups = "groups",
  pair = c("PC1", "PC2"),
  dist_explained = c(45, 22),
  dist_metric = "bray-curtis"
)
```
