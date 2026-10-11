# Color map of a variable

Creates a color map of [setNames](https://rdrr.io/r/stats/setNames.html)
class with names as unique group labels given a column name and values
as hexcode colors. This function is built into the `composition` method
from the abstract class
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) and
inherited by other omics classes, such as;
[metagenomics](https://agusinac.github.io/OmicFlow/reference/metagenomics.md)
and
[proteomics](https://agusinac.github.io/OmicFlow/reference/proteomics.md).

## Usage

``` r
colormap(
  data,
  col_name = lifecycle::deprecated(),
  groups,
  Brewer.palID = "Set2"
)
```

## Arguments

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) or
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html).

- col_name:

  **\[deprecated\]** This argument has been renamed to `groups` for more
  clarity.

- groups:

  A column name of a categorical variable in `data`.

- Brewer.palID:

  A character name that exists in
  [brewer.pal](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  (Default: `"Set2"`).

## Value

A [setNames](https://rdrr.io/r/stats/setNames.html).

## Examples

``` r
library("data.table")
#> 
#> Attaching package: ‘data.table’
#> The following object is masked from ‘package:base’:
#> 
#>     %notin%
dt <- data.table(
  "SAMPLE_ID" = c("sample_1", "sample_2", "sample_3"),
  "treatment" = c("healthy", "tumor", NA)
)

colors <- colormap(data = dt,
                   groups = "treatment")
```
