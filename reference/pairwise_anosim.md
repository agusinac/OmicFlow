# Pairwise ANOSIM

Computes pairwise
[anosim](https://vegandevs.github.io/vegan/reference/anosim.html), given
a distance matrix and a vector of labels. This function is built into
the class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) with
method `ordination()` and inherited by other omics classes, such as;
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md).

## Usage

``` r
pairwise_anosim(
  x,
  groups,
  metadata = NULL,
  perm_design = NULL,
  p.adjust.method = "bonferroni",
  perm = 999,
  ...
)
```

## Arguments

- x:

  A distance matrix in the form of
  [dist](https://rdrr.io/r/stats/dist.html). Obtained from a
  dissimilarity metric, in the case of similarity metric please use
  `1-dist`

- groups:

  A character vector with equal length as the distance matrix.

- metadata:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) as input to the
  function `perm_design` (default: `NULL`).

- perm_design:

  A function that takes the `metadata` and returns a permutation design
  with [how](https://rdrr.io/pkg/permute/man/how.html) (default:
  `NULL`).

- p.adjust.method:

  A character variable to adjust the p-values, see
  [p.adjust](https://rdrr.io/r/stats/p.adjust.html) (default:
  `"bonferroni"`).

- perm:

  A whole number to define the number of permutations in
  [adonis2](https://vegandevs.github.io/vegan/reference/adonis.html)
  (default: `999`).

- ...:

  Additional arguments passed to
  [anosim](https://vegandevs.github.io/vegan/reference/anosim.html).

## Value

A [data.frame](https://rdrr.io/r/base/data.frame.html) containing:

- pairs:

  combinations of group comparisons

- Df:

  Degrees of Freedom

- R2:

  The R squared of the null hypothesis \\H_0\\

- p.value:

  The number of \\F^p\\ higher than the null hypothesis divided by the
  total number of permutations.

- p.adj:

  The adjusted P-value based on the used `p.adjust.method`

## See also

[anosim](https://vegandevs.github.io/vegan/reference/anosim.html)

## Examples

``` r
# Create random data
set.seed(42)
mock_data <- matrix(rnorm(15 * 10), nrow = 15, ncol = 10)

# Create euclidean dissimilarity matrix
mock_dist <- dist(mock_data, method = "euclidean")

# Define group labels, should be equal to number of columns and rows to dist
mock_groups <- rep(c("A", "B", "C"), each = 5)

# Compute pairwise anosim
result <- pairwise_anosim(x = mock_dist, 
                          groups = mock_groups, 
                          p.adjust.method = "bonferroni", 
                          perm = 99)
```
