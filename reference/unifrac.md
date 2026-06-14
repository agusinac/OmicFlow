# Compute UniFrac Dissimilarity from a Dense or Sparse Matrix.

Calculates the UniFrac dissimilarity between samples based on
phylogenetic branch lengths and abundance or presence/absence data.

## Usage

``` r
unifrac(x, tree, weighted = TRUE, normalized = TRUE, threads = 1)
```

## Arguments

- x:

  A [matrix](https://rdrr.io/r/base/matrix.html),
  [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) or
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) of strictly
  positive counts or presence/absence data.

- tree:

  A `phylo` class tree.

- weighted:

  A boolean value, to use abundances (`weighted = TRUE`) or
  absence/presence (`weighted=FALSE`) (default: TRUE).

- normalized:

  A boolean value, whether to normalize weighted UniFrac distances to be
  between 0 and 1 (default: TRUE). Unweighted UniFrac is always
  normalized.

- threads:

  A wholenumber, the number of threads to use in
  [setThreadOptions](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)
  (default: 1).

## Value

A column x column [dist](https://rdrr.io/r/stats/dist.html) object.

## Details

The UniFrac distance between two samples \\A\\ and \\B\\, with
phylogenetic tree edges \\i = 1 \ldots n\\ of lengths \\L_i\\, is
computed differently depending on the `weighted` and `normalized` flags.
When `weighted = FALSE`, input counts are first converted to
presence/absence data.

- Weighted UniFrac (`normalized = FALSE` and `weighted = TRUE`)::

  \\d(A,B) = \frac{\sum\_{i}^n L_i \|A_i - B_i\|}{\sum\_{i}^n L_i (A_i +
  B_i)}\\

- Normalized Weighted UniFrac (`normalized = TRUE` and
  `weighted = TRUE`)::

  \\d(A,B) = \sum\_{i}^n L_i \|A_i - B_i\|\\

- Unweighted UniFrac (`weighted = FALSE`, unweighted is always
  normalized)::

  \\d(A,B) = \frac{\sum\_{i}^n L_i \|A_i - B_i\|}{\sum\_{i}^n L_i
  \max(A_i, B_i)}\\

## References

Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method
for comparing microbial communities. Applied and Environmental
Microbiology, 71(12), 8228–8235.

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

# Weighted UniFrac
unifrac(x = taxa$countData, tree = taxa$treeData, weighted=TRUE, normalized=FALSE)
#>            S100       S103       S115
#> S103 0.38658597                      
#> S115 0.08090148 0.37767607           
#> S120 0.34751952 0.12478228 0.33777195

# Weighted Normalized UniFrac
unifrac(x = taxa$countData, tree = taxa$treeData, weighted=TRUE, normalized=TRUE)
#>           S100      S103      S115
#> S103 0.6314552                    
#> S115 0.1244192 0.6167280          
#> S120 0.5791822 0.2219650 0.5627751

# Unweighted UniFrac
unifrac(x = taxa$countData, tree = taxa$treeData, weighted=FALSE)
#>           S100      S103      S115
#> S103 0.8791970                    
#> S115 0.7630199 0.8686165          
#> S120 0.7981928 0.7444571 0.7713334
```
