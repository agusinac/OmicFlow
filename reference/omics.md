# Superclass omics

OmicFlow is constructed around the main superclass `omics` that in turn
is inherited by
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md),
that offer additional fields or methods. The superclass `omics` contains
both public and private methods, only the public componenents are
documented, whereas the private methods can be accessed via
`omics$private_methods`, which is not supported.

The `omics` class components use
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

## Active bindings

- `metaData`:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  with `SAMPLE_ID` column.

- `featureData`:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  with `FEATURE_ID` column.

- `countData`:

  A [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  with rownames matching `FEATURE_ID` of `featureData` and colnames
  matching `SAMPLE_ID` of `countData`.

## Methods

### Public methods

- [`omics$new()`](#method-omics-initialize)

- [`omics$copy()`](#method-omics-copy)

- [`omics$print()`](#method-omics-print)

- [`omics$reset()`](#method-omics-reset)

- [`omics$removeNAs()`](#method-omics-removeNAs)

- [`omics$feature_subset()`](#method-omics-feature_subset)

- [`omics$sample_subset()`](#method-omics-sample_subset)

- [`omics$samplepair_subset()`](#method-omics-samplepair_subset)

- [`omics$feature_merge()`](#method-omics-feature_merge)

- [`omics$scale()`](#method-omics-scale)

- [`omics$rankstat()`](#method-omics-rankstat)

- [`omics$alpha_diversity()`](#method-omics-alpha_diversity)

- [`omics$composition()`](#method-omics-composition)

- [`omics$distance()`](#method-omics-distance)

- [`omics$ordination()`](#method-omics-ordination)

- [`omics$foldchange()`](#method-omics-foldchange)

- [`omics$autoFlow()`](#method-omics-autoFlow)

- [`omics$clone()`](#method-omics-clone)

------------------------------------------------------------------------

### `omics$new()`

#### Usage

    omics$new(metaData = NULL, countData = NULL, featureData = NULL)

#### Arguments

- `metaData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) (default:
  `NULL`).

- `countData`:

  A path to an existing file,
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html),
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) (default:
  `NULL`).

- `featureData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) (default:
  `NULL`).

#### Returns

A new omics object.

#### Examples

    library("OmicFlow")

    ## Method 1: load from filepath
    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    obj <- omics$new(
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

    obj <- omics$new(
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

    obj <- omics$new(
     metaData = metadata,
     countData = counts
    )

------------------------------------------------------------------------

### `omics$copy()`

This method is very similar to the existing [`$clone()`](#method-clone)
function, except it also creates a new save of the back-up of each class
that is created when [`$new()`](#method-new) is called and can be
reverted to via [`$reset()`](#method-reset)

#### Usage

    omics$copy(deep = FALSE)

#### Arguments

- `deep`:

  A boolean value to create a shallow or deep copy (default: `FALSE`).

#### Returns

A copy of `omics` object.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")

    obj <- omics$new(
     metaData = metadata_file,
     countData = counts_file
    )

    # Perform a modification and copy
    obj$scale()

    cloned <- obj$copy(deep=TRUE)
    cloned$scale(method = "clr")
    cloned$reset() # resets to data after clone creation.

------------------------------------------------------------------------

### `omics$print()`

Displays field parameters of the class via stdout.

#### Usage

    omics$print()

#### Returns

A cli formatted message

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")

    obj <- omics$new(
     metaData = metadata_file,
     countData = counts_file
    )

    # method 1 to call print function
    obj

    # method 2 to call print function
    obj$print()

------------------------------------------------------------------------

### `omics$reset()`

Upon creation of a new class object (i.e. `omics`, `metagenomics` or
`proteomics`) a backup of the original data is created. Since
modification of the object is done by reference and duplicates are not
made, it is possible to revert back to the unmodified data via
`$reset()`. The methods from the super class omics also contains a
private method to prevent any changes to the original object when using
methods, such as `ordination`, `alpha_diversity` or `foldchange`.

#### Usage

    omics$reset()

#### Returns

The object is modified in place and returned (invisibly). If a copy is
required, please use [`$copy()`](#method-copy) beforehand with
`obj2 <- obj$copy()`.

#### Examples

    library(ggplot2)
    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    taxa <- omics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file
    )

    # Performs modifications
    taxa$scale(transform = log2)

    # resets
    taxa$reset()

    # An inbuilt reset function prevents unwanted modification to the taxa object.
    taxa$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))

------------------------------------------------------------------------

### `omics$removeNAs()`

Subset `metaData` via the
[complete.cases](https://rdrr.io/r/stats/complete.cases.html) given a
single or vector of column names, then synchronize the object.

#### Usage

    omics$removeNAs(column)

#### Arguments

- `column`:

  A character vector of a single or multiple column names to remove NAs.

#### Returns

The object is modified in place and returned (invisibly). If a copy is
required, please use [`$copy()`](#method-copy) beforehand with
`obj2 <- obj$copy()`.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    obj$removeNAs(column = "treatment")

------------------------------------------------------------------------

### `omics$feature_subset()`

Subset features using a logical expression evaluated in `featureData`,
then synchronize the object.

#### Usage

    omics$feature_subset(expr)

#### Arguments

- `expr`:

  A logical expression evaluated in the context of `featureData`. Rows
  for which the expression evaluates to `TRUE` are retained.

#### Returns

The object is modified in place and returned (invisibly). If a copy is
required, please use [`$copy()`](#method-copy) beforehand with
`obj2 <- obj$copy()`.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    obj$feature_subset(Genus == "Pseudomonas")

------------------------------------------------------------------------

### `omics$sample_subset()`

Subset samples using a logical expression evaluated in `metaData`, then
synchronize the object.

#### Usage

    omics$sample_subset(expr)

#### Arguments

- `expr`:

  A logical expression evaluated in the context of `metaData`. Rows for
  which the expression evaluates to `TRUE` are retained.

#### Returns

The object is modified in place and returned (invisibly). If a copy is
required, please use [`$copy()`](#method-copy) beforehand with
`obj2 <- obj$copy()`.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    obj$sample_subset(treatment == "tumor")

------------------------------------------------------------------------

### `omics$samplepair_subset()`

Subset samplepairs by choosing the number of unique pairs in `metaData`,
then synchronize the object.

#### Usage

    omics$samplepair_subset(num_unique_pairs = NULL)

#### Arguments

- `num_unique_pairs`:

  An integer value to define the number of pairs to subset. The default
  is NULL, meaning the maximum number of unique pairs will be used to
  subset the data. Let's say you have three samples for each pair, then
  the `num_unique_pairs` will be set to 3.

#### Returns

The object is modified in place and returned (invisibly). If a copy is
required, please use [`$copy()`](#method-copy) beforehand with
`obj2 <- obj$copy()`.

------------------------------------------------------------------------

### `omics$feature_merge()`

Agglomerates features by column, then synchronize the object.

#### Usage

    omics$feature_merge(feature_rank, feature_filter = NULL)

#### Arguments

- `feature_rank`:

  A character (vector) of columns to aggregate from the `featureData`.

- `feature_filter`:

  A character (vector) of characters to remove features via regex
  substring pattern matching (default: `NULL`).

#### Returns

The object is modified in place and returned (invisibly). If a copy is
required, please use [`$copy()`](#method-copy) beforehand with
`obj2 <- obj$copy()`.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    obj$feature_merge(feature_rank = c("Kingdom", "Phylum"))
    obj$feature_merge(feature_rank = "Genus", feature_filter = c("uncultured", "metagenome"))

------------------------------------------------------------------------

### `omics$scale()`

Feature scaling on the `countData`. The `$scale` function is able to
apply transformations element-wise on the positive values, (optional:
add pseudocounts) and perform normalisation or standardisation methods.

#### Usage

    omics$scale(
      method = "tss",
      transform = NULL,
      base = exp(1),
      pseudocount = NULL
    )

#### Arguments

- `method`:

  A character to choose a standardisation/normalisation method, options:
  `"tss"`, `"clr"`, `"binary"`, `"hellinger"`, `"none"` (default:
  `"tss"`).

  `"tss"`

  :   Total-sum scaling. Each column is divided by their column sum,
      converting absolute abundance to relative abundances in range from
      0 to 1.

  `"clr"`

  :   Centered log-ratio-like transformation. Non-zero entries are
      log-transformed using `log` and `base`, then each feature row is
      centered by subtracting by their respective geometric mean. Input
      counts must be strictly positive for his transformation. If
      `pseudocount` is used, then output is equal to
      `vegan::decostand(method = "clr")`, otherwise the current `"clr"`
      approximates the `vegan::decostand(method = "rclr")`, the main
      difference is that vegan replaces all zero's by NAs and then
      excludes them when the arithmetic mean is taken, whereas zero's in
      `"clr"` are only not used for log-transformation but preserved
      when taking the arithmetic mean of each feature row.

  `"binary"`

  :   Presence-absence transformation. Every stored non-zero count is
      replaced by one, while implicit zeros in the sparse matrix remain
      zero.

  `"hellinger"`

  :   Hellinger transformation. Counts are first total-sum scaled with
      `"tss"` to relative abundances, after which the square root of
      each non-zero value is taken.

  `"none"`

  :   No transformation. Typically applied when `transform` is used.

- `transform`:

  A function to apply on the positive values of `countData`, skip
  standardisation/normalisation with `method = "none"` (default:
  `NULL`).

- `base`:

  Input for [log](https://rdrr.io/r/base/Log.html) to use natural
  logarithmic scale, log2, log10 or other (default: `exp(1)`) in CLR.

- `pseudocount`:

  A numeric value to replace zero's (default: `NULL`).

#### Returns

The object is modified in place and returned (invisibly). If a copy is
required, please use [`$copy()`](#method-copy) beforehand with
`obj2 <- obj$copy()`.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )
    # standard relative abundance computation
    obj$scale()

    # CLR
    obj$reset()
    obj$scale(method = "clr")

    # transform
    obj$reset()
    obj$scale(method = "none", transform = log2)

------------------------------------------------------------------------

### `omics$rankstat()`

Rank statistics based on `featureData`. Counts the number of (unique)
features identified for each column of interest from the `featureData`
(excluding the `FEATURE_ID`).

#### Usage

    omics$rankstat(feature_ranks, unique = FALSE)

#### Arguments

- `feature_ranks`:

  A vector of characters that match the `featureData`.

- `unique`:

  A boolean value to display only unique entries in `feature_ranks`
  (default: `FALSE`).

#### Returns

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    plt <- obj$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
    plt

------------------------------------------------------------------------

### `omics$alpha_diversity()`

Alpha diversity computation, statistical testing and visualisation.

It computes the alpha diversity via
[diversity](https://agusinac.github.io/OmicFlow/reference/diversity.md),
performs a pairwise wilcoxon statistical test given `groups` via the
[boxjitter_test](https://agusinac.github.io/OmicFlow/reference/boxjitter_test.md).

#### Usage

    omics$alpha_diversity(
      col_name = lifecycle::deprecated(),
      groups,
      metric = "shannon",
      Brewer.palID = "Set2",
      group_by = lifecycle::deprecated(),
      split_by = NULL,
      evenness = FALSE,
      paired = FALSE,
      p.adjust.method = "fdr"
    )

#### Arguments

- `col_name`:

  **\[deprecated\]** This argument has been renamed to `groups` for more
  clarity.

- `groups`:

  A column name in `metaData` to define the contrast of interest for
  statistical testing and visualisation.

- `metric`:

  An alpha diversity metric as input to
  [diversity](https://agusinac.github.io/OmicFlow/reference/diversity.md),
  options: `"shannon"`, `"simpson"`, `"invsimpson"` (default:
  `"shannon"`).

- `Brewer.palID`:

  A character name for the palette set to be applied, see
  [brewer.pal](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html) or
  [colormap](https://agusinac.github.io/OmicFlow/reference/colormap.md)
  (default: `"Set2"`).

- `group_by`:

  **\[deprecated\]** This argument has been renamed to `split_by` for
  more clarity

- `split_by`:

  A column name to split the groups into chunks for grouped statistical
  test in
  [boxjitter_test](https://agusinac.github.io/OmicFlow/reference/boxjitter_test.md)
  (default: `NULL`).

- `evenness`:

  A boolean wether to divide diversity by number of species, see
  [specnumber](https://vegandevs.github.io/vegan/reference/diversity.html)
  (default: `FALSE`).

- `paired`:

  A boolean value to perform paired analysis in
  [wilcox.test](https://rdrr.io/r/stats/wilcox.test.html) and samplepair
  subsetting via [`$samplepair_subset()`](#method-samplepair_subset)

- `p.adjust.method`:

  A character variable to specify the p.adjust.method to be used
  (default: `"fdr"`).

#### Returns

A list of components:

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) from
  [diversity](https://agusinac.github.io/OmicFlow/reference/diversity.md).

- stats:

  A pairwise statistics from pairwise_wilcox_test.

- plot:

  A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    plt <- obj$alpha_diversity(groups = "treatment",
                               metric = "shannon")

------------------------------------------------------------------------

### `omics$composition()`

Creates a long `data.table` of `SAMPLE_ID`s, the top 10-15 most abundant
bacteria based on the `feature_rank` with any additional columns
specified via `add_cols`. The `$composition()` method works with
aggregated values by running first
[`$feature_merge()`](#method-feature_merge), then sorting for the most
abundant bacteria based on the row sums. The output of `$composition()`
can be used as input to
[composition_plot](https://agusinac.github.io/OmicFlow/reference/composition_plot.md).

#### Usage

    omics$composition(
      feature_rank,
      feature_filter = NULL,
      col_name = lifecycle::deprecated(),
      add_cols = NULL,
      feature_top = 10,
      Brewer.palID = "RdYlBu"
    )

#### Arguments

- `feature_rank`:

  A column name in `featureData` to aggregate via
  [`$feature_merge()`](#method-feature_merge).

- `feature_filter`:

  A character or vector of characters to removes features by regex
  substring pattern matching via
  [`$feature_merge()`](#method-feature_merge) (default: `NULL`).

- `col_name`:

  **\[deprecated\]** This argument has been renamed to `add_cols` for
  more clarity.

- `add_cols`:

  Optional, a column name or vector of multiple column names from
  `metaData` to add to the `data` list name (default": `NULL`).

- `feature_top`:

  A wholenumber of the top features to visualize, the max is 15, due to
  a limit of palettes (default: `10`).

- `Brewer.palID`:

  A character name for the palette set to be applied, see
  [brewer.pal](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html) or
  [colormap](https://agusinac.github.io/OmicFlow/reference/colormap.md)
  (default: `"RdYlBu"`).

#### Returns

A list of components:

- data:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) of
  feature compositions.

- palette:

  A [setNames](https://rdrr.io/r/stats/setNames.html) palette from
  [colormap](https://agusinac.github.io/OmicFlow/reference/colormap.md)
  matching the top features.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    result <- obj$composition(feature_rank = "Genus",
                              feature_filter = c("uncultured"),
                              feature_top = 10)

    plt <- composition_plot(data = result$data,
                            palette = result$palette,
                            feature_rank = "Genus")

------------------------------------------------------------------------

### `omics$distance()`

Compute a distance metric from `countData`

#### Usage

    omics$distance(
      metric,
      weighted = TRUE,
      threads = 1,
      normalize = TRUE,
      base = exp(1)
    )

#### Arguments

- `metric`:

  A dissimilarity metric to be applied on the `countData`, options:

  `"bray"`

  :   [bray](https://agusinac.github.io/OmicFlow/reference/bray.md)

  `"jaccard"`

  :   [jaccard](https://agusinac.github.io/OmicFlow/reference/jaccard.md)

  `"cosine"`

  :   [cosine](https://agusinac.github.io/OmicFlow/reference/cosine.md)

  `"manhattan"`

  :   [manhattan](https://agusinac.github.io/OmicFlow/reference/manhattan.md)

  `"euclidean"`

  :   [euclidean](https://agusinac.github.io/OmicFlow/reference/euclidean.md)

  `"jsd"`

  :   [jsd](https://agusinac.github.io/OmicFlow/reference/jsd.md)

  `"canberra"`

  :   [canberra](https://agusinac.github.io/OmicFlow/reference/canberra.md)

  `"unifrac"`

  :   [unifrac](https://agusinac.github.io/OmicFlow/reference/unifrac.md)

  `"aitchison"`

  :   The `"aitchison"` metric uses the `$scale(method = "clr")` without
      any pseudocounts, see [`$scale()`](#method-scale), followed by
      euclidean dissimilarity computation.

      Ref: Aitchison, J. (1986) The Statistical Analysis of
      Compositional Data. Chapman and Hall, London, 416 p.

- `weighted`:

  A boolean value, to use abundances (`weighted = TRUE`) or
  absence/presence (`weighted=FALSE`) (default: `TRUE`).

- `threads`:

  A wholenumber, indicating the number of threads to use (default: `1`).

- `normalize`:

  A boolean value, whether to normalize weighted UniFrac distances to be
  between 0 and 1. Unweighted UniFrac is always normalized (default:
  `TRUE`).

- `base`:

  Input for [log](https://rdrr.io/r/base/Log.html) to use natural
  logarithmic scale, log2, log10 or other (default: `exp(1)`).

#### Returns

A column x column [dist](https://rdrr.io/r/stats/dist.html) object.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
        metaData = metadata_file,
        countData = counts_file,
        featureData = features_file
    )

    obj$feature_subset(Kingdom == "Bacteria")
    dist <- obj$distance(metric = "bray")

------------------------------------------------------------------------

### `omics$ordination()`

The `$ordination()` method either computes a distance matrix by
specification of the `metric` or takes a custom distance matrix via
`distmat`. It then computes either the pairwise PERMANOVA test on the
defined contrast via `group_by` via
[pairwise_adonis](https://agusinac.github.io/OmicFlow/reference/pairwise_adonis.md)
or the pairwise ANOSIM test via
[pairwise_anosim](https://agusinac.github.io/OmicFlow/reference/pairwise_anosim.md).

Which pairwise test to compute is decided via the `method` argument,
where `"pcoa"` results in the
[pairwise_adonis](https://agusinac.github.io/OmicFlow/reference/pairwise_adonis.md)
and `"nmds"` to
[pairwise_anosim](https://agusinac.github.io/OmicFlow/reference/pairwise_anosim.md).

It also visualises these results as a scree-plot (only when `"pcoa"` is
used), scores-plot (see
[ordination_plot](https://agusinac.github.io/OmicFlow/reference/ordination_plot.md))
and an anova-plot (see
[plot_pairwise_stats](https://agusinac.github.io/OmicFlow/reference/plot_pairwise_stats.md)).
These plots, alongisde the distance matrix, principal components are
returned as a list.

#### Usage

    omics$ordination(
      metric = "bray",
      method = "pcoa",
      group_by,
      distmat = NULL,
      weighted = TRUE,
      normalize = TRUE,
      threads = 1,
      perm_design = NULL,
      perm = 999
    )

#### Arguments

- `metric`:

  A dissimilarity or similarity metric to be applied on the `countData`,
  see [`$distance()`](#method-distance) (default: `"bray"`).

- `method`:

  Ordination method, options: `"pcoa"` or `"nmds"`, see
  [wcmdscale](https://vegandevs.github.io/vegan/reference/wcmdscale.html)
  (default: `"pcoa"`).

  `"pcoa"`

  :   When `"pcoa"` (Principal Coordinate Analysis) is chosen, a
      pairwise PERMANOVA test via
      [pairwise_adonis](https://agusinac.github.io/OmicFlow/reference/pairwise_adonis.md)
      is conducted on the contrast specified by `group_by` with n number
      of permutations, `perm`.

  `"nmds"`

  :   When `"nmds"` (Non-Metric Multidimensional Scaling) is chosen, a
      pairwise ANOSIM test via
      [pairwise_anosim](https://agusinac.github.io/OmicFlow/reference/pairwise_anosim.md)
      is conducted on the contrast specified by `group_by` with n number
      of permutations, `perm`. This method returns 2 dimensions by
      default, and a stress score in the `pcs` list name, that acts as
      the confidence, see
      [metaMDS](https://vegandevs.github.io/vegan/reference/metaMDS.html).
      Since only two dimensions are returned, there is no `scree_plot`
      output.

- `group_by`:

  A character variable in `metaData` to be used for the
  [pairwise_adonis](https://agusinac.github.io/OmicFlow/reference/pairwise_adonis.md)
  or
  [pairwise_anosim](https://agusinac.github.io/OmicFlow/reference/pairwise_anosim.md)
  statistical test.

- `distmat`:

  A custom distance matrix as either
  [dist](https://rdrr.io/r/stats/dist.html) or
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) (default:
  `NULL`).

- `weighted`:

  A boolean value, whether to compute weighted or unweighted
  dissimilarities, see [`$distance()`](#method-distance) (default:
  `TRUE`).

- `normalize`:

  A boolean value, wether to normalize weighted UniFrac distances to be
  between 0 and 1, see
  [unifrac](https://agusinac.github.io/OmicFlow/reference/unifrac.md)
  (default: `TRUE`).

- `threads`:

  A wholenumber, indicating the number of threads to use (default: `1`).

- `perm_design`:

  A function that takes `metaData` and constructs a permutation design
  with [how](https://rdrr.io/pkg/permute/man/how.html), see
  [example](https://agusinac.github.io/OmicFlow/articles/visualisation.html#beta-diversity)
  (default: `NULL`).

- `perm`:

  A whole number to define the number of permutations in
  [adonis2](https://vegandevs.github.io/vegan/reference/adonis.html) and
  [anosim](https://vegandevs.github.io/vegan/reference/anosim.html)
  (default: `999`).

#### Returns

A list of components:

- distmat:

  A distance dissimilarity in
  [matrix](https://rdrr.io/r/base/matrix.html) format.

- stats:

  A statistical test as a
  [data.frame](https://rdrr.io/r/base/data.frame.html).

- pcs:

  principal components as a
  [data.frame](https://rdrr.io/r/base/data.frame.html).

- scree_plot:

  A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object.

- anova_plot:

  A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object.

- scores_plot:

  A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object.

#### Examples

    library("ggplot2")
    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- metagenomics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file,
    )

    pcoa_plots <- obj$ordination(metric = "bray",
                                 method = "pcoa",
                                 group_by = "treatment",
                                 weighted = TRUE)
    pcoa_plots

------------------------------------------------------------------------

### `omics$foldchange()`

Differential feature expression for both paired and non-paired data.

The function performs feature agglomeration via
[`$feature_merge()`](#method-feature_merge) (which can be disabled via
`feature_merge=FALSE`), subsetting to remove NAs in `condition.group`
via [`$removeNAs()`](#method-removeNAs) and finding samplepairs when
`paired = TRUE` is supplied via
[`$subset_samplepairs()`](#method-subset_samplepairs).

The fold-changes can be computed differently based on the `method` and
`aggregate_method` options. Transformations of the data should be done
beforehand via [`$scale()`](#method-scale). Finally, homogeneity of
variance is computed based on the selected `aggregate_method` option. If
`aggregate_method = "mean"` then the
[row_levene](https://rdrr.io/pkg/matrixTests/man/levene.html) is
applied, and if `aggregate_method = "median"` is used then the
[row_brownforsythe](https://rdrr.io/pkg/matrixTests/man/levene.html) is
applied. Any filtering of the results is left to the end-user, volcano
plots are automatically added to the list name `volcano_plot` via
[volcano_plot](https://agusinac.github.io/OmicFlow/reference/volcano_plot.md).

#### Usage

    omics$foldchange(
      condition.group,
      condition_A,
      condition_B,
      method = "identity",
      aggregate_method = "median",
      group_by = lifecycle::deprecated(),
      split_by = NULL,
      feature_merge = FALSE,
      feature_rank = "FEATURE_ID",
      feature_filter = NULL,
      paired = FALSE,
      pvalue.threshold = 0.05,
      logfold.threshold = 0.06,
      abundance.threshold = 0
    )

#### Arguments

- `condition.group`:

  A character of an existing column name in `metaData`, wherein the
  conditions A and B are located.

- `condition_A`:

  A character (vector) in `condition.group`.

- `condition_B`:

  A character (vector) in `condition.group`.

- `method`:

  A character to choose the method of fold-change computation (default:
  `"identity"`).

  `"identity"`

  :   Computes fold-change via `log2(condition_A) - log2(condition_B)`,
      and will handle zero's to prevent `Inf` values. Proportional data
      is also supported and will be automatically detected.

  `"log"`

  :   Computes fold-change via `condition_A - condition_B`.

  How to aggregate the fold-change values for `condition_A` and
  `condition_B` can be controlled with the `aggregate_method` argument.

- `aggregate_method`:

  A function to aggregate the fold-change values from multiple samples
  (default: `"median"`). options:

  `"mean"`

  :   It takes the arithmetic mean for each feature, separately for
      `condition_A` and `condition_B`.

  `"median"`

  :   It takes the median for each feature, separately for `condition_A`
      and `condition_B`.

  `"geomean"`

  :   It computes the geometric mean of the positive counts, similar to
      the estimateSizeFactors, when `type = "poscount"`. The only
      difference, is that the code is adapted to handle a
      [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).

  `"none"`

  :   No aggreation is applied, the whole matrix is used as it is.

- `group_by`:

  **\[deprecated\]** This argument has been renamed to `split_by` for
  more clarity.

- `split_by`:

  A character variable of an existing column in `metaData` to split the
  table in chunks prior to fold-change computation (default: `NULL`).
  When disabled then column names will end with `_in_all`.

- `feature_merge`:

  A boolean value wether to call
  [`$feature_merge()`](#method-feature_merge) (default: `FALSE`).

- `feature_rank`:

  A column in the `featureData` to use as the feature scope (default:
  `"FEATURE_ID"`).

- `feature_filter`:

  A character or vector of characters to remove features via regex
  substring pattern matching via
  [`$feature_merge()`](#method-feature_merge) (default: `NULL`).

- `paired`:

  A boolean value, the paired is only applicable when a `SAMPLEPAIR_ID`
  column exists within the `metaData` (default: `FALSE`), see
  [`$samplepair_subset()`](#method-samplepair_subset).

- `pvalue.threshold`:

  A numeric value used as a p-value threshold to label and color
  significant features (default: `0.05`), see
  [volcano_plot](https://agusinac.github.io/OmicFlow/reference/volcano_plot.md).

- `logfold.threshold`:

  A numeric value used as a fold-change threshold to label and color
  significantly expressed features (default: `0.06`), see
  [volcano_plot](https://agusinac.github.io/OmicFlow/reference/volcano_plot.md).

- `abundance.threshold`:

  A numeric value used as an abundance threshold to size the scatter
  dots based on their mean abundance (default: `0`), see
  [volcano_plot](https://agusinac.github.io/OmicFlow/reference/volcano_plot.md).

#### Returns

A list of components:

- `split_by` subsets:

  A [matrix](https://rdrr.io/r/base/matrix.html) subset from each
  `split_by` separated by `condition_A` and `condition_B`.

- data:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  main output of fold-changes between conditions, contains abundance,
  fold-change, and homogeneity tests.

- volcano_plot:

  A list of
  [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) plots
  for each contrast of `A vs B`, number of plots depend on `split_by`
  and `condition_A` input.

#### Examples

    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- omics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file
    )
    ## Method 1: using non-log transformation
    obj$scale(method = "tss")

    dfe <- obj$foldchange(
     feature_rank = "Genus",
     feature_merge = TRUE,
     method = "identity",
     aggregate_method = "median",
     condition.group = "treatment",
     condition_A = "tumor",
     condition_B = "healthy"
    )

    ## Method 2: using log transformation
    obj$reset()
    obj$scale(method = "clr")

    dfe <- obj$foldchange(
     feature_rank = "Genus",
     feature_merge = TRUE,
     method = "log",
     aggregate_method = "mean",
     condition.group = "treatment",
     condition_A = "tumor",
     condition_B = "healthy"
    )

------------------------------------------------------------------------

### `omics$autoFlow()`

Automated Omics Analysis only works if there are column names in the
`metaData` with the prefix `"CONTRAST_"`, which is detected via
sub-string pattern matching. The transformation of the data is left to
the user via the [`$scale()`](#method-scale) method. The `$autoFlow()`
incorporates the following functions:

- `$rankstat()`:

  This function is only applied when `feature_ranks` is specified,
  useful in
  [metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md),
  see [`$rankstat()`](#method-rankstat) for more information.

- `$alpha_diversity()`:

  This function only works with the `metagenomics` class, uses
  `"shannon"` metric, and works for both paired and non-paired data,
  applied on each `"CONTRAST_"` column, see
  [`$alpha_diversity()`](#method-alpha_diversity) for more information.

- `$ordination()`:

  This function is applied on all classes, multiple `distance_metrics`
  can be supplied, and only the Principal Component method
  (`$ordination(method = "pcoa")`) is used. Additional arguments can be
  controlled with `distance_metrics`, `distmat`, `weighted`, `perm`,
  `threads`, see [`$ordination()`](#method-ordination) for more
  information. **NOTE**: `perm_design` is not supported in
  [`$autoFlow()`](#method-autoFlow).

- `$composition()`:

  This function is applied on the significant contrasts that are
  obtained from either `$alpha_diversity()` or `$ordination()`. Multiple
  `feature_contrasts` can be supplied, by default it only uses the
  `feature_contrasts = "FEATURE_ID"`. The composition plots are created
  by default via the `composition_plot(group_by)` with
  `feature_top = 15`. Additional arguments can be controlled with
  `feature_filter`, see [`$composition()`](#method-composition) for more
  information.

- `$foldchange()`:

  This function is applied on the significant contrasts that are
  obtained from either `$alpha_diversity()` or `$ordination()`. Multiple
  `feature_contrasts` can be supplied, by default it only uses the
  `feature_contrasts = "FEATURE_ID"`. Additional arguments can be
  controlled with `fc_method`, `aggregate_method`, `pvalue.threshold`,
  `logfold.threshold`, `abundance.threshold`, see
  [`$foldchange()`](#method-foldchange) for more information.

The above functions are called in order, only significant contrasts are
further used to create composition and volcano plots. The `autoFlow` is
intended as a first screening of certain contrasts and extra analysis,
such as the use of a permutation design are not supported and should be
applied by the user.

The `report` argument controls if an `rmarkdown` document is created,
and with `filename` the output filename can be specified.

#### Usage

    omics$autoFlow(
      feature_contrast = "FEATURE_ID",
      feature_filter = NULL,
      feature_ranks = NULL,
      distance_metrics = c("bray"),
      distmat = NULL,
      weighted = TRUE,
      fc_method = "log",
      aggregate_method = "median",
      pvalue.threshold = 0.05,
      logfold.threshold = 0.06,
      abundance.threshold = 0,
      perm = 999,
      threads = 1,
      report = TRUE,
      filename = paste0(getwd(), "/report.html")
    )

#### Arguments

- `feature_contrast`:

  A character vector of feature columns in the `featureData` to
  aggregate via [`$feature_merge()`](#method-feature_merge) (default:
  `"FEATURE_ID"`).

- `feature_filter`:

  A character vector to filter unwanted features, see
  [`$feature_merge()`](#method-feature_merge) (default: `NULL`).

- `feature_ranks`:

  A character vector as input to [`$rankstat()`](#method-rankstat)
  (default: `NULL`).

- `distance_metrics`:

  A character vector specifying what (dis)similarity metrics to use in
  [`ordination()`](#method-ordination) (default: `c("bray")`). When you
  are working with log-transformed data it is advised to use the
  `"euclidean"`.

- `distmat`:

  A path to an existing file or a dense/sparse
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) format to be used
  in [`$ordination()`](#method-ordination) (default: `NULL`).

- `weighted`:

  A boolean value, whether to compute weighted or unweighted
  dissimilarities in [`$ordination()`](#method-ordination) (default:
  `TRUE`).

- `fc_method`:

  A character to choose the method in
  [`$foldchange()`](#method-foldchange) (default: `"identity"`).

  `"identity"`

  :   Computes fold-change via `log2(condition_A) - log2(condition_B)`,
      and will handle zero's to prevent `Inf` values. Proportional data
      is also supported and will be automatically detected.

  `"log"`

  :   Computes fold-change via `condition_A - condition_B`.

- `aggregate_method`:

  A function to aggregate the matrix values in
  [`$foldchange()`](#method-foldchange) by taking e.g. the `median` of
  `condition_A` and `condition_B` prior to substraction, or in the case
  of `method = "identity"` to take the median of the fold-change
  `log2(median(A)) - log2(median(B))` (default: `median`).

- `pvalue.threshold`:

  A numeric value used as a p-value threshold to label and color
  significant features in [`$foldchange()`](#method-foldchange)
  (default: `0.05`).

- `logfold.threshold`:

  A numeric value used as a fold-change threshold to label and color
  significantly expressed features in
  [`$foldchange()`](#method-foldchange) (default: `0.06`).

- `abundance.threshold`:

  A numeric value used as an abundance threshold to size the scatter
  dots based on their mean abundance in
  [`$foldchange()`](#method-foldchange) (default: `0`).

- `perm`:

  A wholenumber, number of permutations to compare against the null
  hypothesis of
  [adonis2](https://vegandevs.github.io/vegan/reference/adonis.html) or
  [anosim](https://vegandevs.github.io/vegan/reference/anosim.html) in
  [`$ordination()`](#method-ordination) (default: `999`).

- `threads`:

  Number of threads to use, only used in
  [`$distance()`](#method-distance) when distmat is not supplied
  (default: `1`).

- `report`:

  A boolean value to create a HTML markdown report (default: `FALSE`).
  If `FALSE` a nested list of the plots and data is returned.

- `filename`:

  A character to name the HTML report to be saved in the current working
  directory (default: `paste0(getwd(), "/report.html")`). The
  [`getwd()`](https://rdrr.io/r/base/getwd.html) is required for
  rmarkdown to save it in the right path.

#### Returns

if `report = TRUE` a HTML report is generated, otherwise a list
containing:

- plots:

  - `rankstat_plot` as a
    [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html).

  - `alpha_div_plots` as a list of
    [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html).

  - `pcoa_plots` as a nested matrix of lists with
    [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
    (`contrasts` x `distance_metrics`).

  - `composition_plots` as a nested matrix of lists with
    [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
    (`contrasts` x `feature_contrasts`).

  - `Log2FC_plots` as a nested matrix of lists with
    [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
    (`contrasts` x `feature_contrasts`).

- data:

  - `alpha_div_data` as a list of lists of `data` and `stats`.

  - `pcoa_data` as a nested matrix of lists of `stats`, `dist_mat` and
    `pcs` (`contrasts` x `distance_metrics`).

  - `composition_data` as a nested matrix of lists of `data`
    (`contrasts` x `feature_contrasts`).

  - `Log2FC_data` as a nested matrix of lists of `data` (`contrasts` x
    `feature_contrasts`).

------------------------------------------------------------------------

### `omics$clone()`

The objects of this class are cloneable with this method.

#### Usage

    omics$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `omics$new()`
## ------------------------------------------------

library("OmicFlow")

## Method 1: load from filepath
metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
obj <- omics$new(
 metaData = metadata_file,
 countData = counts_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ countData is loaded.
#> ℹ Created placeholder featureData.

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

obj <- omics$new(
 metaData = metadata,
 featureData = features, # optional
 countData = counts
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ featureData is loaded.
#> ✔ countData is loaded.

## Method 3: you have features and counts in a single data.frame
metadata <- data.frame(SAMPLE_ID = c("S1", "S2", "S3"))
counts <- data.frame(
 S1 = c(2,3,0), 
 S2 = c(2,"",1), 
 S3 = c(2,1,NA), 
 proteins = c("ZEB1", "ZEB2", "MAPK")
)

obj <- omics$new(
 metaData = metadata,
 countData = counts
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ℹ Placeholder featureData extracted from countData.
#> ✔ countData is loaded.

## ------------------------------------------------
## Method `omics$copy()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")

obj <- omics$new(
 metaData = metadata_file,
 countData = counts_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ countData is loaded.
#> ℹ Created placeholder featureData.

# Perform a modification and copy
obj$scale()

cloned <- obj$copy(deep=TRUE)
cloned$scale(method = "clr")
cloned$reset() # resets to data after clone creation.


## ------------------------------------------------
## Method `omics$print()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")

obj <- omics$new(
 metaData = metadata_file,
 countData = counts_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ countData is loaded.
#> ℹ Created placeholder featureData.

# method 1 to call print function
obj
#> 
#> ── <omics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 242 features
#> featureData: 0 attributes × 242 features

# method 2 to call print function
obj$print()
#> 
#> ── <omics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 242 features
#> featureData: 0 attributes × 242 features


## ------------------------------------------------
## Method `omics$reset()`
## ------------------------------------------------

library(ggplot2)
library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

taxa <- omics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ featureData is loaded.
#> ✔ countData is loaded.

# Performs modifications
taxa$scale(transform = log2)

# resets
taxa$reset()

# An inbuilt reset function prevents unwanted modification to the taxa object.
taxa$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))



## ------------------------------------------------
## Method `omics$removeNAs()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

obj$removeNAs(column = "treatment")


## ------------------------------------------------
## Method `omics$feature_subset()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

obj$feature_subset(Genus == "Pseudomonas")
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 3 samples
#> countData: 3 samples × 4 features
#> featureData: 7 attributes × 4 features


## ------------------------------------------------
## Method `omics$sample_subset()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

obj$sample_subset(treatment == "tumor")
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 2 samples
#> countData: 2 samples × 115 features
#> featureData: 7 attributes × 115 features


## ------------------------------------------------
## Method `omics$feature_merge()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

obj$feature_merge(feature_rank = c("Kingdom", "Phylum"))
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 28 features
#> featureData: 7 attributes × 28 features
obj$feature_merge(feature_rank = "Genus", feature_filter = c("uncultured", "metagenome"))
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 20 features
#> featureData: 7 attributes × 20 features


## ------------------------------------------------
## Method `omics$scale()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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
# standard relative abundance computation
obj$scale()

# CLR
obj$reset()
obj$scale(method = "clr")

# transform
obj$reset()
obj$scale(method = "none", transform = log2)


## ------------------------------------------------
## Method `omics$rankstat()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

plt <- obj$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
plt


## ------------------------------------------------
## Method `omics$alpha_diversity()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

plt <- obj$alpha_diversity(groups = "treatment",
                           metric = "shannon")


## ------------------------------------------------
## Method `omics$composition()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

result <- obj$composition(feature_rank = "Genus",
                          feature_filter = c("uncultured"),
                          feature_top = 10)
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 63 features
#> featureData: 7 attributes × 63 features

plt <- composition_plot(data = result$data,
                        palette = result$palette,
                        feature_rank = "Genus")


## ------------------------------------------------
## Method `omics$distance()`
## ------------------------------------------------

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

obj$feature_subset(Kingdom == "Bacteria")
#> 
#> ── <metagenomics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 185 features
#> featureData: 7 attributes × 185 features
dist <- obj$distance(metric = "bray")

## ------------------------------------------------
## Method `omics$ordination()`
## ------------------------------------------------

library("ggplot2")
library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- metagenomics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file,
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

pcoa_plots <- obj$ordination(metric = "bray",
                             method = "pcoa",
                             group_by = "treatment",
                             weighted = TRUE)
#> 'nperm' >= set of all permutations: complete enumeration.
#> Set of permutations < 'minperm'. Generating entire set.
pcoa_plots
#> $dist
#>           S100      S103      S115      S120
#> S100 0.0000000 1.0000000 0.8845188 1.0000000
#> S103 1.0000000 0.0000000 1.0000000 0.9470058
#> S115 0.8845188 1.0000000 0.0000000 1.0000000
#> S120 1.0000000 0.9470058 1.0000000 0.0000000
#> 
#> $anova_data
#>              pairs Df SumsOfSqs   F.Model        R2 p.value p.adj
#> 1 tumor vs healthy  1 0.4197984 0.8395968 0.2956747       1     1
#> 
#> $pcs
#>           PC1           PC2           PC3 treatment samples
#>         <num>         <num>         <num>    <char>  <char>
#> 1:  0.3808548  0.000000e+00  4.422594e-01     tumor       1
#> 2: -0.3808548 -4.735029e-01 -3.124747e-16     tumor       2
#> 3:  0.3808548  1.937010e-16 -4.422594e-01   healthy       3
#> 4: -0.3808548  4.735029e-01 -1.430586e-16   healthy       4
#> 
#> $scree_plot

#> 
#> $anova_plot

#> 
#> $scores_plot
#> Too few points to calculate an ellipse
#> Too few points to calculate an ellipse
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_path()`).

#> 


## ------------------------------------------------
## Method `omics$foldchange()`
## ------------------------------------------------

library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

obj <- omics$new(
 metaData = metadata_file,
 countData = counts_file,
 featureData = features_file
)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
#> ✔ featureData is loaded.
#> ✔ countData is loaded.
## Method 1: using non-log transformation
obj$scale(method = "tss")

dfe <- obj$foldchange(
 feature_rank = "Genus",
 feature_merge = TRUE,
 method = "identity",
 aggregate_method = "median",
 condition.group = "treatment",
 condition_A = "tumor",
 condition_B = "healthy"
)
#> 
#> ── <omics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 65 features
#> featureData: 7 attributes × 65 features
#> Warning: row_wilcoxon_twosample: 54 of the rows had ties: cannot compute exact p-values with ties.
#> First occurrence at row 3
#> Warning: row_brownforsythe: 65 of the rows had no groups with at least 3 observations.
#> First occurrence at row 1

## Method 2: using log transformation
obj$reset()
obj$scale(method = "clr")

dfe <- obj$foldchange(
 feature_rank = "Genus",
 feature_merge = TRUE,
 method = "log",
 aggregate_method = "mean",
 condition.group = "treatment",
 condition_A = "tumor",
 condition_B = "healthy"
)
#> 
#> ── <omics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 65 features
#> featureData: 7 attributes × 65 features
#> Warning: row_wilcoxon_twosample: 54 of the rows had ties: cannot compute exact p-values with ties.
#> First occurrence at row 3
#> Warning: row_levene: 65 of the rows had no groups with at least 3 observations.
#> First occurrence at row 1
```
