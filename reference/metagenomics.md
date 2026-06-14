# Sub-class metagenomics

This is a sub-class that is compatible to data obtained from either 16S
rRNA marker-gene sequencing or shot-gun metagenomics sequencing. It
inherits all methods from the abstract class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) and only
adapts the `initialize` function. It supports BIOM format data (v2.1.0
from <http://biom-format.org/>) in both HDF5 and JSON format, also
pre-existing data structures can be used or text files. When omics data
is very large, data loading becomes very expensive. It is therefore
recommended to use the [`reset()`](#method-reset) method to reset your
changes. Every omics class creates an internal memory efficient back-up
of the data, the resetting of changes is an instant process.

## See also

[omics](https://agusinac.github.io/OmicFlow/reference/omics.md)

[volcano_plot](https://agusinac.github.io/OmicFlow/reference/volcano_plot.md)

## Super class

[`omics`](https://agusinac.github.io/OmicFlow/reference/omics.md) -\>
`metagenomics`

## Active bindings

- `treeData`:

  A "phylo" class, see
  [as.phylo](https://rdrr.io/pkg/ape/man/as.phylo.html).

## Methods

### Public methods

- [`metagenomics$new()`](#method-metagenomics-initialize)

- [`metagenomics$write_biom()`](#method-metagenomics-write_biom)

- [`metagenomics$foldchange()`](#method-metagenomics-foldchange)

- [`metagenomics$clone()`](#method-metagenomics-clone)

Inherited methods

- [`omics$alpha_diversity()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-alpha_diversity)
- [`omics$autoFlow()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-autoFlow)
- [`omics$composition()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-composition)
- [`omics$copy()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-copy)
- [`omics$distance()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-distance)
- [`omics$feature_merge()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-feature_merge)
- [`omics$feature_subset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-feature_subset)
- [`omics$ordination()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-ordination)
- [`omics$print()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-print)
- [`omics$rankstat()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-rankstat)
- [`omics$removeNAs()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-removeNAs)
- [`omics$reset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-reset)
- [`omics$sample_subset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-sample_subset)
- [`omics$samplepair_subset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-samplepair_subset)
- [`omics$scale()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-scale)
- [`omics$validate()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-validate)

------------------------------------------------------------------------

### `metagenomics$new()`

Initializes the metagenomics class object with `metagenomics$new()`

#### Usage

    metagenomics$new(
      countData = NULL,
      metaData = NULL,
      featureData = NULL,
      treeData = NULL,
      biomData = NULL,
      feature_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    )

#### Arguments

- `countData`:

  A path to an existing file or a dense/sparse
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) format.

- `metaData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  data.frame.

- `featureData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  data.frame.

- `treeData`:

  A path to an existing newick file or class "phylo", see
  [read.tree](https://rdrr.io/pkg/ape/man/read.tree.html).

- `biomData`:

  A path to an existing biom file, version 2.1.0
  (http://biom-format.org/), see
  [h5read](https://huber-group-embl.github.io/rhdf5/reference/h5_read.html).

- `feature_names`:

  A character vector to name the feature names that fit the supplied
  `featureData`.

#### Returns

A new `metagenomics` object.

------------------------------------------------------------------------

### `metagenomics$write_biom()`

Creates a BIOM file in HDF5 format of the loaded items via
['new()'](#method-new), which is compatible to the python biom-format
version 2.1, see http://biom-format.org.

#### Usage

    metagenomics$write_biom(filename)

#### Arguments

- `filename`:

  A character variable of either the full path of filename of the biom
  file (e.g. `output.biom`)

#### Examples

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

    taxa$write_biom(filename = "output.biom")
    file.remove("output.biom")

------------------------------------------------------------------------

### `metagenomics$foldchange()`

Differential feature expression (DFE) Total Sample Sum (TSS) transformed
values for both paired and non-paired test.

The function performs feature agglomeration, subsetting to remove NAs in
`condition.group`, finding samplepairs and finally feature scaling prior
to fold-change computation. Based on the `transform` method,
fold-changes will be computed either via subtraction or division.

#### Usage

    metagenomics$foldchange(
      condition.group,
      condition_A,
      condition_B,
      group_by = NULL,
      feature_rank = "FEATURE_ID",
      feature_filter = NULL,
      paired = FALSE,
      normalize = FALSE,
      pvalue.threshold = 0.05,
      logfold.threshold = 0.06,
      abundance.threshold = 0
    )

#### Arguments

- `condition.group`:

  A character variable of an existing column name in `metaData`, wherein
  the conditions A and B are located.

- `condition_A`:

  A character value or vector of characters.

- `condition_B`:

  A character value or vector of characters.

- `group_by`:

  A character variable of an existing column in `metaData` to split the
  table in chunks prior to fold-change computation (default: NULL).

- `feature_rank`:

  A character or vector of characters in the `featureData` to aggregate
  via [`feature_merge()`](#method-feature_merge) (default:
  `"FEATURE_ID"`).

- `feature_filter`:

  A character or vector of characters to remove features via regex
  pattern (default: `NULL`).

- `paired`:

  A boolean value, the paired is only applicable when a `SAMPLEPAIR_ID`
  column exists within the `metaData`. See
  [wilcox.test](https://rdrr.io/r/stats/wilcox.test.html) and
  [`samplepair_subset()`](#method-samplepair_subset).

- `normalize`:

  A boolean value wether to normalize via Total Sample Sums (TSS) or not
  (default: `FALSE`).

- `pvalue.threshold`:

  A numeric value used as a p-value threshold to label and color
  significant features (default: 0.05).

- `logfold.threshold`:

  A numeric value used as a fold-change threshold to label and color
  significantly expressed features (default: 0.06).

- `abundance.threshold`:

  A numeric value used as an abundance threshold to size the scatter
  dots based on their mean abundance (default: 0.01).

#### Returns

- `data` A long
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  table.

- `volcano_plot` A
  [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

- `A` A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  table for (each) condition A.

- `B` A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  table for (each) condition B.

#### Examples

    library("ggplot2")
    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file
    )

    dfe <- obj$foldchange(feature_rank = "Genus",
                          paired = FALSE,
                          condition.group = "treatment",
                          condition_A = c("healthy"),
                          condition_B = c("tumor"))

------------------------------------------------------------------------

### `metagenomics$clone()`

The objects of this class are cloneable with this method.

#### Usage

    metagenomics$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `metagenomics$write_biom()`
## ------------------------------------------------

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

taxa$write_biom(filename = "output.biom")
file.remove("output.biom")
#> [1] TRUE


## ------------------------------------------------
## Method `metagenomics$foldchange()`
## ------------------------------------------------

library("ggplot2")
library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ featureData is loaded.
#> ✔ countData is loaded.
#> ℹ Final steps .. cleaning & creating back-up
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 242 features
#> featureData: 7 attributes × 242 features

dfe <- obj$foldchange(feature_rank = "Genus",
                      paired = FALSE,
                      condition.group = "treatment",
                      condition_A = c("healthy"),
                      condition_B = c("tumor"))
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 64 features
#> featureData: 7 attributes × 64 features
```
