# Compute Jensen-Shannon Divergence from a Dense or Sparse Matrix.

Calculates the Jensen-Shannon divergence of a Matrix pairwise for each
column.

## Usage

``` r
jsd(x, weighted = TRUE, threads = 1)
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

The Jensen-Shannon divergence between two probability distributions
\\A\\ and \\B\\, each of length \\n\\, is defined as:

\\ d(A, B) = \frac{1}{2} D\_{KL}(A \parallel M) + \frac{1}{2} D\_{KL}(B
\parallel M) \\

where \\M = \frac{1}{2} (A + B)\\ is the mixture distribution, and
\\D\_{KL}\\ is the Kullback-Leibler divergence. When weighted is set to
FALSE, counts are replaced by presence/absence data.

## References

Lin, J. (1991). Divergence measures based on the Shannon entropy. IEEE
Transactions on Information Theory, 37(1), 145-151.

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

jsd(taxa$countData)
#>           S100      S103      S115
#> S103 1.0000000                    
#> S115 0.8135729 1.0000000          
#> S120 1.0000000 0.9196102 1.0000000
```
