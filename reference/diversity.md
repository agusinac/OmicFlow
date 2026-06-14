# Sparse implementation of Alpha Diversity Metrics

Computes the alpha diversity based on Shannon index, simpson or
invsimpson. Code is adapted from
[diversity](https://vegandevs.github.io/vegan/reference/diversity.html)
and uses
[sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) in
triplet format over the dense matrix. The code is much faster and memory
efficient, while still being mathematically correct. This function is
built into the class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) with
method `alpha_diversity()` and inherited by other omics classes, such
as;
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md).

## Usage

``` r
diversity(
  x,
  metric = c("shannon", "simpson", "invsimpson"),
  normalize = TRUE,
  base = exp(1)
)
```

## Arguments

- x:

  A [matrix](https://rdrr.io/r/base/matrix.html),
  [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) or
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html).

- metric:

  A character variable for metric; shannon, simpson or invsimpson.

- normalize:

  A boolean variable for sample normalization by column sums.

- base:

  Input for [log](https://rdrr.io/r/base/Log.html) to use natural
  logarithmic scale, log2, log10 or other (default: `exp(1)`).

## Value

A numeric vector with type double.

## See also

[diversity](https://vegandevs.github.io/vegan/reference/diversity.html)

## Examples

``` r
n_row <- 1000
n_col <- 100
density <- 0.2
num_entries <- n_row * n_col
num_nonzero <- round(num_entries * density)

set.seed(123)
positions <- sample(num_entries, num_nonzero, replace=FALSE)
row_idx <- ((positions - 1) %% n_row) + 1
col_idx <- ((positions - 1) %/% n_row) + 1

values <- runif(num_nonzero, min = 0, max = 1)
sparse_mat <- sparseMatrix(
  i = row_idx,
  j = col_idx,
  x = values,
  dims = c(n_row, n_col)
)

# Alpha diversity is computed on column level
## Transpose the sparseMatrix if required with t() from Matrix R package.
result <- OmicFlow::diversity(
  x = sparse_mat,
  metric = "shannon"
) 
```
