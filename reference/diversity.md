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
diversity(x, metric = "shannon", normalize = TRUE, base = exp(1))
```

## Arguments

- x:

  A [matrix](https://rdrr.io/r/base/matrix.html),
  [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) or
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html).

- metric:

  A character to specify what metric to use, options: `"shannon"`,
  `"simpson"`, `"invsimpson"` (default: `"shannon"`).

  `"shannon"`

  :   Shannon Diversity Index, \\H = -\sum_i p_i
      \log\_{\code{base}}(p_i)\\, where \\p_i\\ is the abundance or
      relative abundance of feature \\i\\. When `normalize = TRUE`,
      counts are converted to within-sample relative abundances before
      calculation. The logarithm base is controlled by `base`; the
      natural logarithm is commonly used.

  `"simpson"`

  :   Simpson Diversity Index, \\1 - \sum_i p_i^2\\. Higher values
      indicate greater diversity, with values near zero indicating
      domination by one feature. When `normalize = TRUE`, \\p_i\\
      denotes within-sample relative abundance.

  `"invsimpson"`

  :   Inverse Simpson Index, \\1 / \sum_i p_i^2\\. When calculated from
      relative abundances, this is the effective number of equally
      abundant features that would yield the observed concentration.
      Values increase with richness and evenness.

- normalize:

  A boolean variable to normalise the each value by their column sum
  (default: `TRUE`).

- base:

  A numeric variable in [log](https://rdrr.io/r/base/Log.html) to use
  natural logarithmic scale, log2, log10 or other in `"shannon"`
  (default: `exp(1)`).

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
sparse_mat <- Matrix::sparseMatrix(
  i = row_idx,
  j = col_idx,
  x = values,
  dims = c(n_row, n_col)
)

# Alpha diversity is computed on column level
## Transpose the sparseMatrix if required with `t()` from Matrix R package.
result <- OmicFlow::diversity(
  x = sparse_mat,
  metric = "shannon"
) 
```
