# Compute Manhattan Dissimilarity from a Dense or Sparse Matrix.

Calculates the Manhattan dissimilarity of a Matrix pairwise for each
column.

## Usage

``` r
manhattan(x, weighted = TRUE, threads = 1)
```

## Arguments

- x:

  A [matrix](https://rdrr.io/r/base/matrix.html),
  [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) or
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html).

- weighted:

  A boolean value, to use abundances (`weighted = TRUE`) or
  absence/presence (`weighted=FALSE`) (default: TRUE).

- threads:

  A wholenumber, the number of threads to use in
  [setThreadOptions](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)
  (default: 1).

## Value

A column x column [dist](https://rdrr.io/r/stats/dist.html) object.

## Details

The Manhattan dissimilarity between two samples \\A\\ and \\B\\, each of
length \\n\\, is defined as:

\\d(A, B) = \sum\_{i}^n \|A_i - B_i\|\\

where \\A_i\\ and \\B_i\\ are the abundances of the \\i\\-th feature in
sample \\A\\ and \\B\\, respectively. When weighted is set to FALSE,
counts are replaced by presence/absence data.

## References

Deza, M. M., & Deza, E. (2009). Encyclopedia of Distances. Springer
Science & Business Media., 313.

## Examples

``` r
library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")

taxa <- metagenomics$new(
    metaData = metadata_file,
    countData = counts_file,
    featureData = features_file,
    treeData = tree_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ featureData is loaded.
#> ✔ countData is loaded.
#> ✔ treeData is loaded.
#> ℹ Final steps .. cleaning & creating back-up
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 242 features
#> featureData: 7 attributes × 242 features
#> treeData: 242 tips × 241 nodes

taxa$feature_subset(Kingdom == "Bacteria")
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 185 features
#> featureData: 7 attributes × 185 features
#> treeData: 185 tips × 184 nodes
taxa$scale(method = "tss")

manhattan(taxa$countData)
#>          S100     S103     S115
#> S103 2.000000                  
#> S115 1.733627 2.000000         
#> S120 2.000000 1.845190 2.000000
```
