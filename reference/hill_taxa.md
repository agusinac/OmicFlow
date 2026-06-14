# Sparse implementation of Hill numbers

Computes the hill numbers for q is 0, 1 or 2. Code is adapted from
hill_taxa and uses
[sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) in
triplet format over the dense matrix. The code is much faster and memory
efficient, while still being mathematical correct.

## Usage

``` r
hill_taxa(x, q = 0, normalize = TRUE, base = exp(1))
```

## Arguments

- x:

  A [matrix](https://rdrr.io/r/base/matrix.html),
  [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) or
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html).

- q:

  A wholenumber for 0, 1 or 2, default is 0.

- normalize:

  A boolean variable for sample normalization by column sums.

- base:

  Input for [log](https://rdrr.io/r/base/Log.html) to use natural
  logarithmic scale, log2, log10 or other.

## Value

A numeric vector with type double.

## See also

hill_taxa

## Examples

``` r
library("Matrix")

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

result <- OmicFlow::hill_taxa(
 x = sparse_mat,
 q = 2
)
```
