# Diversity plot

Creates an Alpha diversity plot. This function is built into the class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) with
method `alpha_diversity()`. It computes the pairwise wilcox test, paired
or non-paired, given a data frame and adds useful labelling.

## Usage

``` r
diversity_plot(
  data,
  values,
  col_name,
  group_by = NULL,
  palette,
  method,
  paired = FALSE,
  p.adjust.method = "fdr"
)
```

## Arguments

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) or
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  computed from
  [diversity](https://agusinac.github.io/OmicFlow/reference/diversity.md).

- values:

  A column name of a continuous variable.

- col_name:

  A column name of a categorical variable.

- group_by:

  A column name to perform grouped statistical test (default: NULL).

- palette:

  An object with names and hexcode or color names, see
  [colormap](https://agusinac.github.io/OmicFlow/reference/colormap.md).

- method:

  A character variable indicating what method is used to compute the
  diversity.

- paired:

  A boolean value to perform paired analysis in
  [wilcox.test](https://rdrr.io/r/stats/wilcox.test.html).

- p.adjust.method:

  A character variable to specify the p.adjust.method to be used
  (Default: fdr).

## Value

A
[ggplot2](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object to be further modified

## Examples

``` r
library("ggplot2")
 
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

div <- OmicFlow::diversity(
  x = sparse_mat,
  metric = "shannon"
)

dt <- data.table::data.table(
  "shannon" = div,
  "treatment" = c(rep("healthy", n_col / 2), rep("tumor", n_col / 2)),
  "sex" = c(rep("male", n_col / 4), rep("female", n_col / 4))
)

colors <- OmicFlow::colormap(dt, "treatment")

# Comparing two groups
plt <- OmicFlow::diversity_plot(
 data = dt,
 values = "shannon",
 col_name = "treatment",
 palette = colors,
 method = "shannon",
 paired = FALSE,
 p.adjust.method = "fdr"
)

# Performing a test while stratifying the plot in two groups
plt <- OmicFlow::diversity_plot(
 data = dt,
 values = "shannon",
 col_name = "treatment",
 group_by = "sex",
 palette = colors,
 method = "shannon",
 paired = FALSE,
 p.adjust.method = "fdr"
)
```
