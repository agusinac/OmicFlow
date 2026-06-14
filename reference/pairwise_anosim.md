# Pairwise anosim (ANOSIM) computation

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
  perm = 999
)
```

## Arguments

- x:

  A distance matrix in the form of
  [dist](https://rdrr.io/r/stats/dist.html). Obtained from a
  dissimilarity metric, in the case of similarity metric please use
  `1-dist`

- groups:

  A vector (column from a table) of labels.

- metadata:

  A data.table or data.frame of extra metadata for `perm_design`
  (default: NULL).

- perm_design:

  A function that takes a data.frame and constructs a permutation design
  with [how](https://rdrr.io/pkg/permute/man/how.html) (default: NULL).

- p.adjust.method:

  P adjust method see [p.adjust](https://rdrr.io/r/stats/p.adjust.html)

- perm:

  Number of permutations to compare against the null hypothesis of
  anosim (default: `perm=999`).

## Value

A [data.frame](https://rdrr.io/r/base/data.frame.html) of

- pairs that are used

- R2 of H_0

- p value of F^p \> F

- p adjusted

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
