# Compute Jaccard Dissimilarity from a Dense or Sparse Matrix.

Calculates the Jaccard dissimilarity of a Matrix pairwise for each
column.

## Usage

``` r
jaccard(x, weighted = TRUE, threads = 1)
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

The weighted Jaccard disimilarity between two samples \\A\\ and \\B\\,
each of length \\n\\, is defined as:

\\d(A,B) = 1 - \frac{ \sum\_{i}^{n} \min(A_i, B_i) }{ \sum\_{i}^{n}
\max(A_i, B_i) }\\

where \\A_i\\ and \\B_i\\ are the abundances of the \\i\\-th feature in
sample \\A\\ and \\B\\, respectively. When weighted is set to FALSE,
abundances are changed to 1 (classical Jaccard for binary data).

## References

Jaccard, P. (1912) The distribution of the flora in the alpine zone. New
Phytologist, 11(2), 37–50. library("OmicFlow")

metadata_file \<- system.file("extdata", "metadata.tsv", package =
"OmicFlow") counts_file \<- system.file("extdata", "counts.tsv", package
= "OmicFlow") features_file \<- system.file("extdata", "features.tsv",
package = "OmicFlow") tree_file \<- system.file("extdata",
"tree.newick", package = "OmicFlow")

taxa \<- metagenomics\$new( metaData = metadata_file, countData =
counts_file, featureData = features_file, treeData = tree_file )

taxa\$feature_subset(Kingdom == "Bacteria") taxa\$scale(method = "tss")

jaccard(taxa\$countData)
