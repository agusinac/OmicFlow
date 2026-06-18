# Abstract omics class

This is the abstract class 'omics', contains a variety of methods that
are inherited and applied in the omics classes:
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md).

## Details

Every class is created with the
[R6Class](https://r6.r-lib.org/reference/R6Class.html) method. Methods
are either public or private, and only the public components are
inherited by other omic classes. The omics class by default uses a
[sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) and
[data.table](https://rdrr.io/pkg/data.table/man/data.table.html) data
structures for quick and efficient data manipulation and returns the
object by reference, same as the R6 class. The method by reference is
very efficient when dealing with big data.

## References

Aitchison, J. (1986) The Statistical Analysis of Compositional Data.
Chapman and Hall, London, 416 p.

## See also

[diversity_plot](https://agusinac.github.io/OmicFlow/reference/diversity_plot.md)

[composition_plot](https://agusinac.github.io/OmicFlow/reference/composition_plot.md)

[bray](https://agusinac.github.io/OmicFlow/reference/bray.md),
[canberra](https://agusinac.github.io/OmicFlow/reference/canberra.md),
[cosine](https://agusinac.github.io/OmicFlow/reference/cosine.md),
[jaccard](https://agusinac.github.io/OmicFlow/reference/jaccard.md),
[jsd](https://agusinac.github.io/OmicFlow/reference/jsd.md),
[manhattan](https://agusinac.github.io/OmicFlow/reference/manhattan.md),
[unifrac](https://agusinac.github.io/OmicFlow/reference/unifrac.md)

[ordination_plot](https://agusinac.github.io/OmicFlow/reference/ordination_plot.md),
[plot_pairwise_stats](https://agusinac.github.io/OmicFlow/reference/plot_pairwise_stats.md),
[pairwise_anosim](https://agusinac.github.io/OmicFlow/reference/pairwise_anosim.md),
[pairwise_adonis](https://agusinac.github.io/OmicFlow/reference/pairwise_adonis.md)

[volcano_plot](https://agusinac.github.io/OmicFlow/reference/volcano_plot.md)

## Active bindings

- `metaData`:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  with `SAMPLE_ID` column.

- `featureData`:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  with `FEATURE_ID` column.

- `countData`:

  A dense or sparse
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html).

## Methods

### Public methods

- [`omics$new()`](#method-omics-initialize)

- [`omics$copy()`](#method-omics-copy)

- [`omics$validate()`](#method-omics-validate)

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

Wrapper function that is inherited and adapted for each omics class. The
omics classes requires a metadata samplesheet, that is validated by the
metadata_schema.json. It requires a column `SAMPLE_ID` and optionally a
`SAMPLEPAIR_ID` can be supplied. The `SAMPLE_ID` will be used to link
the metaData to the countData, and will act as the key during subsetting
of other columns. To create a new object use
[[`new()`](https://rdrr.io/r/methods/new.html)](#method-new) method. Do
notice that the abstract class only checks if the metadata is valid! The
`countData` and `featureData` will not be checked, these are handled by
the sub-classes. Using the omics class to load your data is not
supported and still experimental.

#### Usage

    omics$new(countData = NULL, featureData = NULL, metaData = NULL)

#### Arguments

- `countData`:

  A path to an existing file or a dense/sparse
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) format.

- `featureData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  data.frame.

- `metaData`:

  A path to an existing file,
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  data.frame.

#### Returns

A new `omics` object.

------------------------------------------------------------------------

### `omics$copy()`

Create a copy of the object-class

This method is very similar to the existing [`clone()`](#method-clone)
function, except it also resets the back-up of the OmicFlow data types
that is invoked with [`reset()`](#method-reset)

#### Usage

    omics$copy(deep = FALSE)

#### Arguments

- `deep`:

  A boolean value to create a shallow or deep copy.

#### Returns

A copy of `omics` object

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

### `omics$validate()`

Validates an input metadata against the JSON schema. The metadata should
look as follows and should not contain any empty spaces. For example;
`'sample 1'` is not allowed, whereas `'sample1'` is allowed!

Acceptable column headers:

- SAMPLE_ID (required)

- SAMPLEPAIR_ID (optional)

- CONTRAST\_ (optional), used for [`autoFlow()`](#method-autoFlow).

- VARIABLE\_ (optional), not supported yet.

This function is used during the creation of a new object via
[[`new()`](https://rdrr.io/r/methods/new.html)](#method-new) to validate
the supplied metadata via a filepath or existing
[data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
[data.frame](https://rdrr.io/r/base/data.frame.html).

#### Usage

    omics$validate()

#### Returns

None

------------------------------------------------------------------------

### `omics$print()`

Displays parameters of the omics class via stdout.

#### Usage

    omics$print()

#### Returns

object in place

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

Upon creation of a new `omics` object a small backup of the original
data is created. Since modification of the object is done by reference
and duplicates are not made, it is possible to `reset` changes to the
class. The methods from the abstract class omics also contains a private
method to prevent any changes to the original object when using methods
such as `ordination` `alpha_diversity` or `foldchange`.

#### Usage

    omics$reset()

#### Returns

object in place

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

Remove NAs from `metaData` and updates the `countData`.

#### Usage

    omics$removeNAs(column)

#### Arguments

- `column`:

  The column from where NAs should be removed, this can be either a
  wholenumbers or characters. Vectors are also supported.

#### Returns

object in place

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

Feature subset (based on `featureData`), automatically applies data
synchronization.

#### Usage

    omics$feature_subset(...)

#### Arguments

- `...`:

  Expressions that return a logical value, and are defined in terms of
  the variables in `featureData`. Only rows for which all conditions
  evaluate to TRUE are kept.

#### Returns

object in place

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

Sample subset (based on `metaData`), automatically applies
synchronization.

#### Usage

    omics$sample_subset(...)

#### Arguments

- `...`:

  Expressions that return a logical value, and are defined in terms of
  the variables in `metaData`. Only rows for which all conditions
  evaluate to TRUE are kept.

#### Returns

object in place

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

Samplepair subset (based on `metaData`), automatically applies
synchronization.

#### Usage

    omics$samplepair_subset(num_unique_pairs = NULL)

#### Arguments

- `num_unique_pairs`:

  An integer value to define the number of pairs to subset. The default
  is NULL, meaning the maximum number of unique pairs will be used to
  subset the data. Let's say you have three samples for each pair, then
  the `num_unique_pairs` will be set to 3.

#### Returns

object in place

------------------------------------------------------------------------

### `omics$feature_merge()`

Agglomerates features by column, automatically applies synchronization.

#### Usage

    omics$feature_merge(feature_rank, feature_filter = NULL)

#### Arguments

- `feature_rank`:

  A character value or vector of columns to aggregate from the
  `featureData`.

- `feature_filter`:

  A character value or vector of characters to remove features via regex
  pattern.

#### Returns

object in place

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

Feature scaling on the `countData`. The `scale` function is able to
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
  `tss`, `clr`, `binary`, `hellinger`, `none` (default: `"tss"`).

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

object in place

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

Rank statistics based on `featureData`

#### Usage

    omics$rankstat(feature_ranks, unique = FALSE)

#### Arguments

- `feature_ranks`:

  A vector of characters or integers that match the `featureData`.

- `unique`:

  A boolean value to display only unique entries in `feature_ranks`.

#### Details

Counts the number of features identified for each column, for example in
case of 16S metagenomics it would be the number of OTUs or ASVs on
different taxonomy levels.

#### Returns

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

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

    plt <- obj$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
    plt

------------------------------------------------------------------------

### `omics$alpha_diversity()`

Alpha diversity based on
[diversity](https://agusinac.github.io/OmicFlow/reference/diversity.md)

#### Usage

    omics$alpha_diversity(
      col_name,
      metric = c("shannon", "invsimpson", "simpson"),
      Brewer.palID = "Set2",
      group_by = NULL,
      evenness = FALSE,
      paired = FALSE,
      p.adjust.method = "fdr"
    )

#### Arguments

- `col_name`:

  A character variable from the `metaData`.

- `metric`:

  An alpha diversity metric as input to
  [diversity](https://agusinac.github.io/OmicFlow/reference/diversity.md).

- `Brewer.palID`:

  A character name for the palette set to be applied, see
  [brewer.pal](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html) or
  [colormap](https://agusinac.github.io/OmicFlow/reference/colormap.md).

- `group_by`:

  A column name to perform grouped statistical test in
  [diversity_plot](https://agusinac.github.io/OmicFlow/reference/diversity_plot.md)
  (default: NULL).

- `evenness`:

  A boolean wether to divide diversity by number of species, see
  [specnumber](https://vegandevs.github.io/vegan/reference/diversity.html).

- `paired`:

  A boolean value to perform paired analysis in
  [wilcox.test](https://rdrr.io/r/stats/wilcox.test.html) and samplepair
  subsetting via [`samplepair_subset()`](#method-samplepair_subset)

- `p.adjust.method`:

  A character variable to specify the p.adjust.method to be used,
  default is 'fdr'.

#### Returns

A list of components:

- `div` A [data.frame](https://rdrr.io/r/base/data.frame.html) from
  [diversity](https://agusinac.github.io/OmicFlow/reference/diversity.md).

- `stats` A pairwise statistics from
  [pairwise_wilcox_test](https://rpkgs.datanovia.com/rstatix/reference/wilcox_test.html).

- `plot` A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
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

    plt <- obj$alpha_diversity(col_name = "treatment",
                               metric = "shannon")

------------------------------------------------------------------------

### `omics$composition()`

Creates a table most abundant compositional features. Also assigns a
color blind friendly palette for visualizations.

#### Usage

    omics$composition(
      feature_rank,
      feature_filter = NULL,
      col_name = NULL,
      feature_top = c(10, 15),
      Brewer.palID = "RdYlBu"
    )

#### Arguments

- `feature_rank`:

  A character variable in `featureData` to aggregate via
  [`feature_merge()`](#method-feature_merge).

- `feature_filter`:

  A character or vector of characters to removes features by regex
  pattern.

- `col_name`:

  Optional, a character or vector of characters to add to the final
  compositional data output.

- `feature_top`:

  A wholenumber of the top features to visualize, the max is 15, due to
  a limit of palettes.

- `Brewer.palID`:

  A character name for the palette set to be applied, see
  [brewer.pal](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html) or
  [colormap](https://agusinac.github.io/OmicFlow/reference/colormap.md).

#### Returns

A list of components:

- `data` A
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) of
  feature compositions.

- `palette` A [setNames](https://rdrr.io/r/stats/setNames.html) palette
  from
  [colormap](https://agusinac.github.io/OmicFlow/reference/colormap.md).

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
      normalized = TRUE,
      base = exp(1)
    )

#### Arguments

- `metric`:

  A dissimilarity metric to be applied on the `countData`, thus far
  supports 'bray', 'jaccard', 'cosine', 'manhattan', 'aitchison',
  'euclidean', 'jsd' (jensen-shannon divergence), 'canberra' and
  'unifrac' when a tree is provided via `treeData`, see
  [`distance()`](#method-distance).

- `weighted`:

  A boolean value, to use abundances (`weighted = TRUE`) or
  absence/presence (`weighted=FALSE`) (default: TRUE).

- `threads`:

  A wholenumber, indicating the number of threads to use (Default: 1).

- `normalized`:

  A boolean value, whether to normalize weighted UniFrac distances to be
  between 0 and 1. Unweighted UniFrac is always normalized (default:
  TRUE).

- `base`:

  Input for [log](https://rdrr.io/r/base/Log.html) to use natural
  logarithmic scale, log2, log10 or other (default: `exp(1)`).

- `pseudocount`:

  A numeric value to replace zero's, used in
  [[`scale()`](https://rdrr.io/r/base/scale.html)](#method-scale)
  (default: `1e-15`).

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

Ordination of `countData` with statistical testing.

#### Usage

    omics$ordination(
      metric = "bray",
      method = c("pcoa", "nmds"),
      group_by,
      distmat = NULL,
      weighted = TRUE,
      threads = 1,
      perm_design = NULL,
      perm = 999
    )

#### Arguments

- `metric`:

  A dissimilarity or similarity metric to be applied on the `countData`,
  thus far supports 'bray', 'jaccard', 'cosine', 'manhattan', 'jsd'
  (jensen-shannon divergence), 'canberra' and 'unifrac' when a tree is
  provided via `treeData`, see [`distance()`](#method-distance).

- `method`:

  Ordination method, supports "pcoa" and "nmds", see
  [wcmdscale](https://vegandevs.github.io/vegan/reference/wcmdscale.html).

- `group_by`:

  A character variable in `metaData` to be used for the
  [pairwise_adonis](https://agusinac.github.io/OmicFlow/reference/pairwise_adonis.md)
  or
  [pairwise_anosim](https://agusinac.github.io/OmicFlow/reference/pairwise_anosim.md)
  statistical test.

- `distmat`:

  A custom distance matrix in either
  [dist](https://rdrr.io/r/stats/dist.html) or
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) format.

- `weighted`:

  A boolean value, whether to compute weighted or unweighted
  dissimilarities (default: `TRUE`).

- `threads`:

  A wholenumber, indicating the number of threads to use (Default: 1).

- `perm_design`:

  A function that takes `metaData` and constructs a permutation design
  with [how](https://rdrr.io/pkg/permute/man/how.html) (default:
  `NULL`).

- `perm`:

  A wholenumber, number of permutations to compare against the null
  hypothesis of
  [adonis2](https://vegandevs.github.io/vegan/reference/adonis.html) and
  [anosim](https://vegandevs.github.io/vegan/reference/anosim.html)
  (default: `perm=999`).

#### Returns

A list of components:

- `distmat` A distance dissimilarity in
  [matrix](https://rdrr.io/r/base/matrix.html) format.

- `stats` A statistical test as a
  [data.frame](https://rdrr.io/r/base/data.frame.html).

- `pcs` principal components as a
  [data.frame](https://rdrr.io/r/base/data.frame.html).

- `scree_plot` A
  [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

- `anova_plot` A
  [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

- `scores_plot` A
  [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

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

Differential feature expression (DFE) on log-transformed values for both
paired and non-paired test.

The function performs feature agglomeration, subsetting to remove NAs in
`condition.group` and finding samplepairs. It expects that the data is
already log-transformed, this can be accomplished via
[[`scale()`](https://rdrr.io/r/base/scale.html)](#method-scale)

#### Usage

    omics$foldchange(
      condition.group,
      condition_A,
      condition_B,
      group_by = NULL,
      feature_rank = "FEATURE_ID",
      feature_filter = NULL,
      paired = FALSE,
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
  table in chunks prior to fold-change computation (default: `NULL`).
  When disabled then column names will end with `_in_all`.

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
  table for (each) condition A

- `B` A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  table for (each) condition B

#### Examples

    library("ggplot2")
    library("OmicFlow")

    metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

    obj <- omics$new(
     metaData = metadata_file,
     countData = counts_file,
     featureData = features_file
    )
    obj$scale(method = "clr")

    dfe <- obj$foldchange(feature_rank = "Genus",
                          paired = FALSE,
                          condition.group = "treatment",
                          condition_A = c("healthy"),
                          condition_B = c("tumor"))

------------------------------------------------------------------------

### `omics$autoFlow()`

Automated Omics Analysis based on the `metaData`, see
[`validate()`](#method-validate). For now only works with headers that
start with prefix `CONTRAST_`. If the data is from the class `omics` or
`proteomics` than FDR adjusted p-values are computed for the volcano
plots. Log-transformed values will lead to the skipping of
[`composition()`](#method-composition) and
[`alpha_diversity()`](#method-alpha_diversity) methods.

#### Usage

    omics$autoFlow(
      feature_contrast = "FEATURE_ID",
      feature_filter = NULL,
      feature_ranks = NULL,
      distance_metrics = c("bray"),
      distmat = NULL,
      weighted = TRUE,
      pvalue.threshold = 0.05,
      logfold.threshold = 1,
      abundance.threshold = 0.01,
      perm = 999,
      threads = 1,
      report = TRUE,
      filename = paste0(getwd(), "/report.html")
    )

#### Arguments

- `feature_contrast`:

  A character vector of feature columns in the `featureData` to
  aggregate via [`feature_merge()`](#method-feature_merge) (default:
  `"FEATURE_ID"`).

- `feature_filter`:

  A character vector to filter unwanted features, (default: `NULL`).

- `feature_ranks`:

  A character vector as input to [`rankstat()`](#method-rankstat)
  (default: `NULL`).

- `distance_metrics`:

  A character vector specifying what (dis)similarity metrics to use
  (default: `c("bray")`) When you are working with log-transformed data
  it is advised to use the `euclidean`.

- `distmat`:

  A path to an existing file or a dense/sparse
  [Matrix](https://rdrr.io/pkg/Matrix/man/Matrix.html) format (default:
  `NULL`).

- `weighted`:

  A boolean value, whether to compute weighted or unweighted
  dissimilarities (default: `TRUE`).

- `pvalue.threshold`:

  A numeric value, the p-value is used to include/exclude composition
  and foldchanges plots coming from alpha- and beta diversity analysis
  (default: 0.05).

- `logfold.threshold`:

  A numeric value used as a fold-change threshold to label and color
  significantly expressed features, see
  [`foldchange()`](#method-foldchange) (Default: 1).

- `abundance.threshold`:

  A numeric value used as an abundance threshold to size the scatter
  dots based on their mean abundance, see
  [`foldchange()`](#method-foldchange) (default: 0.01).

- `perm`:

  A wholenumber, number of permutations to compare against the null
  hypothesis of
  [adonis2](https://vegandevs.github.io/vegan/reference/adonis.html) or
  [anosim](https://vegandevs.github.io/vegan/reference/anosim.html)
  (default: 999).

- `threads`:

  Number of threads to use, only used in
  [`distance()`](#method-distance) when distmat is not supplied
  (default: 1).

- `report`:

  A boolean value to create a HTML markdown report (default: `FALSE`).
  If `FALSE` a nested list of the plots and data is returned.

- `filename`:

  A character to name the HTML report to be saved in the current working
  directory (default: `paste0(getwd(), "/report.html")`). The
  [`getwd()`](https://rdrr.io/r/base/getwd.html) is required for
  rmarkdown to save it in the right path.

#### Returns

List of plots/data or rendered HTML report

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
#> ! Created placeholder featureData.

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
#> ! Created placeholder featureData.

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

plt <- obj$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
plt


## ------------------------------------------------
## Method `omics$alpha_diversity()`
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

plt <- obj$alpha_diversity(col_name = "treatment",
                           metric = "shannon")


## ------------------------------------------------
## Method `omics$composition()`
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
#>           PC1           PC2           PC3  groups samples
#>         <num>         <num>         <num>  <char>  <char>
#> 1:  0.3808548  0.000000e+00  4.422594e-01   tumor       1
#> 2: -0.3808548 -4.735029e-01 -3.124747e-16   tumor       2
#> 3:  0.3808548  1.937010e-16 -4.422594e-01 healthy       3
#> 4: -0.3808548  4.735029e-01 -1.430586e-16 healthy       4
#> 
#> $scree_plot
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_col()`).

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

library("ggplot2")
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
obj$scale(method = "clr")

dfe <- obj$foldchange(feature_rank = "Genus",
                      paired = FALSE,
                      condition.group = "treatment",
                      condition_A = c("healthy"),
                      condition_B = c("tumor"))
#> 
#> ── <omics> object 
#> metaData: 9 variables × 4 samples
#> countData: 4 samples × 64 features
#> featureData: 7 attributes × 64 features
```
