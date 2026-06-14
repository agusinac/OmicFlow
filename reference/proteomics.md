# Sub-class proteomics

This is a sub-class that is compatible to preprocessed data obtained
from https://fragpipe.nesvilab.org/. It inherits all methods from the
abstract class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) and only
adapts the `initialize` function. It supports pre-existing data
structures or paths to text files. When omics data is very large, data
loading becomes very expensive. It is therefore recommended to use the
[`reset()`](#method-reset) method to reset your changes. Every omics
class creates an internal memory efficient back-up of the data, the
resetting of changes is an instant process.

## See also

[omics](https://agusinac.github.io/OmicFlow/reference/omics.md)

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

- [`proteomics$clone()`](#method-proteomics-clone)

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
- [`omics$validate()`](https://agusinac.github.io/OmicFlow/reference/omics.html#method-validate)

------------------------------------------------------------------------

### `proteomics$new()`

Initializes the proteomics class object with `proteomics$new()`

#### Usage

    proteomics$new(
      countData = NULL,
      metaData = NULL,
      featureData = NULL,
      treeData = NULL
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

#### Returns

A new `proteomics` object.

------------------------------------------------------------------------

### `proteomics$clone()`

The objects of this class are cloneable with this method.

#### Usage

    proteomics$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
