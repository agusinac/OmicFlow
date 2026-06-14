# Color map of a variable

Creates an object of hexcode colors with names given a vector of
characters. This function is built into the `ordination` method from the
abstract class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) and
inherited by other omics classes, such as;
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md).

## Usage

``` r
colormap(data, col_name, Brewer.palID = "Set2")
```

## Arguments

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) or
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html).

- col_name:

  A column name of a categorical variable.

- Brewer.palID:

  A character name that exists in
  [brewer.pal](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  (Default: `"Set2"`).

## Value

A [setNames](https://rdrr.io/r/stats/setNames.html).

## Examples

``` r
library("data.table")
dt <- data.table(
  "SAMPLE_ID" = c("sample_1", "sample_2", "sample_3"),
  "treatment" = c("healthy", "tumor", NA)
)

colors <- colormap(data = dt,
                   col_name = "treatment")
```
