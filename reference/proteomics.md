# Sub-class proteomics

OmicFlow is constructed around the main superclass
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) that in
turn is inherited by
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and proteomics, that offer additional fields or methods. The superclass
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) contains
both public and private methods, only the public componenents are
documented, whereas the private methods can be accessed via
`proteomics$private_methods`, which is not supported.

The [omics](https://agusinac.github.io/OmicFlow/reference/omics.md)
class components use
[Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) and
[data.table](https://rdrr.io/pkg/data.table/man/data.table.html) in the
background for fast loading and data wrangling, these data structures
return the data by reference, so does the R6 class, therefore classes in
OmicFlow do not create copies of the object unless called via `$copy()`.

All classes needs to be first initialised via the `$new()` method and
require the `metaData` and `countData` field components. Alternatively,
you can also initialise with only the `metaData` and later add the other
fields via the active binding, but this will not create a back-up and is
not advised! For more hands-on examples see [Getting Started with
OmicFlow](https://agusinac.github.io/OmicFlow/articles/getting-started.html).

### Additional arguments in `proteomics$new()`

The `proteomics` class also offers the loading of a phylogenetic tree as
a `phylo` class (see
[as.phylo](https://rdrr.io/pkg/ape/man/as.phylo.html)) that can be
supplied to `treeData`. When a `treeData` is supplied than every class
field is arranged according to the tree tip labels.

## Input requirements

Every class in OmicFlow requires at least a `metaData` and `countData`
field to be specified. The `metaData` must contain a `SAMPLE_ID` and
optionally a `SAMPLEPAIR_ID`, see section "Metadata validation" below.

The `countData` can be both a matrix or data frame. If the original
contains rownames, then these will become the `FEATURE_ID` in
`featureData` (automatically created), otherwise a random feature
identifier is assigned. If the `countData` is a filetype or a data frame
with both matching sample identifiers in `SAMPLE_ID` and other
additional columns, then these will be automatically split, columns
matching the `SAMPLE_ID` will be placed in the `countData`, any other
columns will be placed in front of the `FEATURE_ID` column in
`featureData`.

## Metadata validation

Every class in OmicFlow begins by validating the sample metadata. This
validation checks that sample identifiers match the sample IDs present
in the abundance tables.

The validation rules are defined in JSON format in the omics class.
Metadata can be supplied as either:

- A CSV or TSV file.

- A `data.table`.

In both cases, metadata must include a header row and follow a
one-row-per-sample structure. For file inputs, the header must be on the
first line. Columns not described below are permitted and are ignored
during metadata validation.

### Minimum requirement

The metadata must contain a `SAMPLE_ID` column.

- Each row must have a unique, non-empty sample identifier.

- Sample identifiers cannot contain spaces; use underscores (`_`) or
  dashes (`-`) instead.

Example metadata:


    SAMPLE_ID  SAMPLEPAIR_ID  CONTRAST_Treatment  VARIABLE_Age
    S1         P1             Drug                42
    S2         P1             Placebo             36
    S3         P2             Drug                51

### Column types and naming rules

**Required column**

|  |  |  |
|----|----|----|
| Column | Type | Rules |
| `SAMPLE_ID` | string | Unique, non-empty, no spaces; one identifier per sample row |

**Optional standard columns**

|  |  |  |
|----|----|----|
| Column | Type | Rules |
| `SAMPLEPAIR_ID` | string | Optional; no spaces. Used for samples that are paired or originate from the same individual or subject. |

Some functions require a `SAMEPLEPAIR_ID` in the `metaData` when
`paired = TRUE` is used. The user can add or remove the `SAMPLEPAIR_ID`
at anytime to the `metaData` via active binding.

**Pattern-based columns**

Additional variables can be defined with special column-name prefixes:

- `CONTRAST_...`: Grouping or category labels for differential
  comparisons. For example, `CONTRAST_Treatment` may contain `Drug` and
  `Placebo`.

- `VARIABLE_...`: Numeric or string variables for statistical analysis.
  For example, `VARIABLE_Age` may contain values such as `42` and `51`.

Pattern-based columns are only used in
[`$autoFlow()`](#method-autoFlow). Currently, only columns with the
`CONTRAST_` prefix are supported.

## Super class

[`omics`](https://agusinac.github.io/OmicFlow/reference/omics.md) -\>
`proteomics`

## Active bindings

- `treeData`:

  A "phylo" class, see
  [as.phylo](https://rdrr.io/pkg/ape/man/as.phylo.html).

## Methods

### Public methods

- [`proteomics$new()`](#method-proteomics-initialize)

Inherited methods

- [`omics$alpha_diversity()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-alpha_diversity)
- [`omics$autoFlow()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-autoFlow)
- [`omics$composition()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-composition)
- [`omics$copy()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-copy)
- [`omics$distance()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-distance)
- [`omics$feature_merge()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-feature_merge)
- [`omics$feature_subset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-feature_subset)
- [`omics$foldchange()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-foldchange)
- [`omics$ordination()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-ordination)
- [`omics$print()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-print)
- [`omics$rankstat()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-rankstat)
- [`omics$removeNAs()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-removeNAs)
- [`omics$reset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-reset)
- [`omics$sample_subset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-sample_subset)
- [`omics$samplepair_subset()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-samplepair_subset)
- [`omics$scale()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-scale)

------------------------------------------------------------------------

### `proteomics$new()`

#### Usage

    proteomics$new(
      countData = NULL,
      metaData = NULL,
      featureData = NULL,
      treeData = NULL
    )

#### Arguments

- `countData`:

  A path to an existing file,
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html),
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) (default:
  `NULL`).

- `metaData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) (default:
  `NULL`).

- `featureData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) (default:
  `NULL`).

- `treeData`:

  A path to an existing newick file or class "phylo", see
  [read.tree](https://rdrr.io/pkg/ape/man/read.tree.html) (default:
  `NULL`).

#### Returns

A new proteomics object.

#### Examples

    library("OmicFlow")

    ## Method 1: load from filepath
    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    obj <- proteomics$new(
     metaData = metadata_file,
     countData = counts_file
    )

    ## Method 2: Load from data.frame or matrix
    n_cols <- 5
    n_rows <- 100
    n_vals <- n_cols * n_rows
    metadata <- data.frame("SAMPLE_ID" = paste0("Sample_", 1:n_cols))
    features <- data.frame("FEATURE_ID" = paste0("protein_", 1:n_rows))
    counts <- Matrix::Matrix(
     1:n_vals, nrow = n_rows, ncol = n_cols,
     dimnames = list(features$FEATURE_ID, metadata$SAMPLE_ID)
    )

    obj <- proteomics$new(
     metaData = metadata,
     featureData = features, # optional
     countData = counts
    )

    ## Method 3: you have features and counts in a single data.frame
    metadata <- data.frame(SAMPLE_ID = c("S1", "S2", "S3"))
    counts <- data.frame(
     S1 = c(2,3,0),
     S2 = c(2,"",1),
     S3 = c(2,1,NA),
     proteins = c("ZEB1", "ZEB2", "MAPK")
    )

    obj <- proteomics$new(
     metaData = metadata,
     countData = counts
    )

## Examples

``` r

## ------------------------------------------------
## Method `proteomics$new()`
## ------------------------------------------------

library("OmicFlow")

## Method 1: load from filepath
metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
obj <- proteomics$new(
 metaData = metadata_file,
 countData = counts_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ countData is loaded.
#> ℹ Created placeholder featureData.
#> 
#> ── <proteomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 242 features
#> featureData: 0 attributes × 242 features

## Method 2: Load from data.frame or matrix
n_cols <- 5
n_rows <- 100
n_vals <- n_cols * n_rows
metadata <- data.frame("SAMPLE_ID" = paste0("Sample_", 1:n_cols))
features <- data.frame("FEATURE_ID" = paste0("protein_", 1:n_rows))
counts <- Matrix::Matrix(
 1:n_vals, nrow = n_rows, ncol = n_cols, 
 dimnames = list(features$FEATURE_ID, metadata$SAMPLE_ID)
)

obj <- proteomics$new(
 metaData = metadata,
 featureData = features, # optional
 countData = counts
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ featureData is loaded.
#> ✔ countData is loaded.
#> 
#> ── <proteomics> object 
#> metaData: 1 variables × 5 samples
#> countData: 5 samples × 100 features
#> featureData: 0 attributes × 100 features

## Method 3: you have features and counts in a single data.frame
metadata <- data.frame(SAMPLE_ID = c("S1", "S2", "S3"))
counts <- data.frame(
 S1 = c(2,3,0), 
 S2 = c(2,"",1), 
 S3 = c(2,1,NA), 
 proteins = c("ZEB1", "ZEB2", "MAPK")
)

obj <- proteomics$new(
 metaData = metadata,
 countData = counts
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ℹ Placeholder featureData extracted from countData.
#> ✔ countData is loaded.
#> 
#> ── <proteomics> object 
#> metaData: 1 variables × 3 samples
#> countData: 3 samples × 3 features
#> featureData: 1 attributes × 3 features
```
