# Sub-class metagenomics

OmicFlow is constructed around the main superclass
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) that in
turn is inherited by metagenomics and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md),
that offer additional fields or methods. The superclass
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) contains
both public and private methods, only the public componenents are
documented, whereas the private methods can be accessed via
`metagenomics$private_methods`, which is not supported.

The [omics](https://agusinac.github.io/OmicFlow/reference/omics.md)
class components use
[Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) and
[data.table](https://rdrr.io/pkg/data.table/man/data.table.html) in the
background for fast loading and data wrangling, these data structures
return the data by reference, so does the R6 class, therefore classes in
OmicFlow do not create copies of the object unless called via `$copy()`.

All classes needs to be first initialised via the `$new()` method and
require the `metaData` and `countData` (or `biomData` in case of
`metagenomics`) field components. Alternatively, you can also initialise
with only the `metaData` and later add the other fields via the active
binding, but this will not create a back-up and is not advised! For more
hands-on examples see [Getting Started with
OmicFlow](https://agusinac.github.io/OmicFlow/articles/getting-started.html).

### Additional arguments in `metagenomics$new()`

The `metagenomics` class also offers the loading of a phylogenetic tree
as a `phylo` class (see
[as.phylo](https://rdrr.io/pkg/ape/man/as.phylo.html)) that can be
supplied to `treeData`. When a `treeData` is supplied than every class
field is arranged according to the tree tip labels.

The `biomData` argument expects a [BIOM
v2.1.0.](http://biom-format.org/) in either HDF5 or JSON format. The
`biomData` must be supplied via a filepath and required contents are
checked upfront, such as the existence of
` "/observation/matrix/data", "/observation/matrix/indptr", "/observation/matrix/indices", "/observation/metadata/taxonomy", "/observation/ids", "/sample/ids" `.

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
`metagenomics`

## Active bindings

- `treeData`:

  A "phylo" class, see
  [as.phylo](https://rdrr.io/pkg/ape/man/as.phylo.html).

## Methods

### Public methods

- [`metagenomics$new()`](#method-metagenomics-initialize)

- [`metagenomics$write_biom()`](#method-metagenomics-write_biom)

- [`metagenomics$clone()`](#method-metagenomics-clone)

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

### `metagenomics$new()`

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

- `biomData`:

  A path to an existing biom file, [version
  2.1.0](http://biom-format.org/) (see
  [h5read](https://huber-group-embl.github.io/rhdf5/reference/h5_read.html))
  or JSON format (see
  [validate_json_file](https://coolbutuseless.github.io/package/yyjsonr/reference/validate_json_file.html))
  (default: `NULL`).

- `feature_names`:

  A character vector to name the feature names that fit the supplied
  `featureData` (default:
  `c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")`).

#### Returns

A new metagenomics object.

#### Examples

    library("OmicFlow")

    ## Method 1: Load from filepath
    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")

    taxa <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,  # optional
     treeData = tree_file          # optional
    )

    ## Method 2: Load from a biom file
    biom_file <- system.file("extdata", "biomv2.biom", package = "OmicFlow")

    taxa <- metagenomics$new(
     metaData = metadata_file,
     biomData = biom_file
    )

------------------------------------------------------------------------

### `metagenomics$write_biom()`

Creates a BIOM file in HDF5 format, which is compatible to the [python
biom-format v2.1.0.](http://biom-format.org) It only saves the
`featureData` and `countData` to the biom format, no `treeData` is
supported thus far.

#### Usage

    metagenomics$write_biom(filename)

#### Arguments

- `filename`:

  A character variable of a non-existing file path (e.g. `"output.biom"`
  or `"my/path/output.biom"`)

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
## Method `metagenomics$new()`
## ------------------------------------------------

library("OmicFlow")

## Method 1: Load from filepath
metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")

taxa <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,  # optional
 treeData = tree_file          # optional
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

## Method 2: Load from a biom file
biom_file <- system.file("extdata", "biomv2.biom", package = "OmicFlow")

taxa <- metagenomics$new(
 metaData = metadata_file,
 biomData = biom_file
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
```
