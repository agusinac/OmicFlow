# OmicFlow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.4.0 - [2025-10-23]

### `Added`
- [#14](https://github.com/agusinac/OmicFlow/issues/14) `unifrac` for both weighted/unweighted and normalized options.
- [#14](https://github.com/agusinac/OmicFlow/issues/14) `bray` for both weighted/unweighted options.
- [#14](https://github.com/agusinac/OmicFlow/issues/14) `jaccard` for both weighted/unweighted options.
- [#14](https://github.com/agusinac/OmicFlow/issues/14) `jsd` for both weighted/unweighted options.
- [#14](https://github.com/agusinac/OmicFlow/issues/14) `cosine` for both weighted/unweighted options.
- [#14](https://github.com/agusinac/OmicFlow/issues/14) `canberra` for both weighted/unweighted options.
- [#14](https://github.com/agusinac/OmicFlow/issues/14) `manhattan` for both weighted/unweighted options.
- `omics$distance` as wrapper for all new dissimilarity metrics.
- [#18](https://github.com/agusinac/OmicFlow/issues/18) `pairwise_adonis2` and `pairwise_anosim` now allow custom permutation designs from `permute` R package.
- Unit tests for dissimilarity metrics

### `Changed`
- `metagenomics` and `proteomics` classes now align the countData (rows) and featureData by the treeData tip.labels.
- `omics$ordination` uses now `omics$distance` to compute dissimilarity metrics (previously computed via `rbiom`).
- Example data loaded from `inst/extdata` instead of re-using `mock.rds` object.

### `Fixed`
- [#15](https://github.com/agusinac/OmicFlow/issues/15) Fixed flipped log2 values when condition A or B is zero.
- [#17](https://github.com/agusinac/OmicFlow/issues/17) Docker is built via `pak` from github instead of CRAN.

### `Deprecated`
- [#13](https://github.com/agusinac/OmicFlow/issues/13) Removed `viridis` dependency, `purrr:map` is replaced by `lapply`
- Removed `rbiom` and it's dependency `slam`.

## v1.3.2 - [2025-09-20]

### `Deprecated`
- `testthat` removed autoflow test

## v1.3.1 - [2025-09-03]

### `Added`
- rhub yaml to check for additional platforms

### `Fixed`
- `testthat` now uses seed & `tempdir()` to meet CRAN policy
- `report.md` now also displays volcano plots
- `omics$autoFlow()` saves report html in correct current path

## v1.3.0 - [2025-08-13]

### `Added`
- Docker support for R package `autoflow` [#10](https://github.com/agusinac/OmicFlow/issues/10)

### `Changed`
- `autoFlow` and `report.md` now contain downloadable data [#11](https://github.com/agusinac/OmicFlow/issues/11)
- Renamed `feature_glom` to `feature_merge`
- Added `@returns` documentation in manuals
- Updated `citations.md`
- Finalized `readme.md`

### `Fixed`
- working directory is set to current working directory in `rmarkdown::render` (compatible in docker/bioconda etc.)

## v1.2.1 - [2025-08-05]

### `Added`
- Included github actions workflow for `covr` code coverage and R CMD build check.
- github MIT license
- created folder `inst/extdata` with test data

### `Changed`
- testthat now only checks presence of `write_biom` and `autoFlow` files.
- All other testthat functions are replaced by snapshot testing, issues [#2](https://github.com/agusinac/OmicFlow/issues/2) and [#3](https://github.com/agusinac/OmicFlow/issues/3)

### `Fixed`
- `on.exit` now solves issue during error occurance that class items are changed, issue [#3](https://github.com/agusinac/OmicFlow/issues/3).

## v1.2.0 - [2025-08-03]

### `Added`
- included `exec` folder with `autoFlow.R` function to be called from the command line [#8](https://github.com/agusinac/OmicFlow/issues/8).
- New `private` functions `check_matrix` and `check_table` for loading of files or existing data structures, issue [#6](https://github.com/agusinac/OmicFlow/issues/6).

### `Changed`
- `autoFlow.R` now has option to select sub-class from the command line.

## v1.1.0 - [2025-08-02]

### `Added`
- included `inst` folder with `report.Rmd` and `css.styles` 
- `autoFlow` now automatically loads `report.Rmd` from `inst` folder [#3](https://github.com/agusinac/OmicFlow/issues/3).
- Test displayed on plots changes based on boolean paired value.

### `Changed`
- refactoring of `autoflow`, handles now both paired and unpaired data [#3](https://github.com/agusinac/OmicFlow/issues/3).
- `combine_conditions` automatically check if comparison is significant and updates the list with new conditions.
- updated manuals
- `foldchange.R` now handles zero's and doesn't return Inf values

## v1.0.7 - [2025-07-29]

### `Added`
- `volcano_plot` now offers option to visualize only abundant bacteria

### `Fixed`
- Wrong order of normalization and feature aggregation [#7](https://github.com/agusinac/OmicFlow/issues/7)

## v1.0.6 - [2025-07-28]

### `Changed`
- Updated documentation according to issues [#2](https://github.com/agusinac/OmicFlow/issues/2) and [#4](https://github.com/agusinac/OmicFlow/issues/4)

### `Fixed`
- Improved regex pattern in json file, now checks and doesnt allow for spaces.

### `Deprecated`
- removed personal files, cleaned up directory.

## v1.0.5 - [2025-07-27]

### `Added`
- Included `@examples` in documentation [#4](https://github.com/agusinac/OmicFlow/issues/4).

## v1.0.4 - [2025-07-24]

### `Fixed`
- Performed unit-testing [#2](https://github.com/agusinac/OmicFlow/issues/2), [#3](https://github.com/agusinac/OmicFlow/issues/3).

## v1.0.3 - [2025-07-23]

### `Added`
- Created a `citation.md` for all used R packages in bibtex format
- Added more error handling and documentation to functions/methods [#2](https://github.com/agusinac/OmicFlow/issues/2).

### `Changed`
- Updated documentation, manual [#2](https://github.com/agusinac/OmicFlow/issues/2)

### `Deprecated`
- Removed non-used functions

## v1.0.2 - [2025-07-22]

### `Added`
- created `column_exist` for efficient error handling of missing/empty columns in tables

### `Fixed`
- `write_biom` is now functional and tested compared to python API `biom-format`.
- `feature_filter` is not compatible with newer version of `feature_glom`
- `featureData` now replaces empty strings with NA, compatible with `metaData`.

## v1.0.1 - [2025-07-18]

### `Fixed`
- `write_biom` is now functional and tested compared to python API `biom-format`

## v1.0.0 - [2025-07-17]

### `Added`
- `samplepair_subset` finds automatically pairs and subsets class.
- created metadata validation based on json schema, function `validate`.
- additional error handling is added with cli R package.

### `Changed`
- `metataxonomics` is now called `metagenomics` sub-class, supports both hdf5, json formats.
- `metagenomics` class first loads `omics` class and then fills up any additional information.
- ordination plot now only shows ecclipses on T distribution.

### `Deprecated`
- Removed `find_pairs`

## v0.9.5 - [2025-06-11]

### `Changed`
- `unpaired_fold` and `paired_fold` are now combined into a single function `foldchange`.

### `Deprecated`
- Removed `doParallel` and `foreach` from foldchange computation.

## v0.9.4 - [2025-06-03]

### `Changed`
- `feature_glom` now handles multiple columns or single ones.
- `metataxonomics` sub-class now automatically renames the last taxonomy columns.

## v0.9.3 - [2025-05-08]

### `Added`
- `read_sparseTable` now efficiently reads in tsv, txt, csv, url or compressed files cleans it before creating a `sparseMatrix`

## v0.9.2 - [2025-03-27]

### `Added`
- error-handling for treeData

### `Changed`
- Extra check for metadata & countdata alignment based of sample-ids
- Improved visualizations of triplot, using ggplot functions
- `tools` class is now called `omics` class.
- `diversity_plot` now only shows significant pvalues, supports p.adjust methods

## v0.9.1 - [2025-03-11]

### `Changed`
- `feature_glom` can now be repeated multiple times without throwing an error.
- Replaced `ggpubr::compare_means` by `rstatix::pairwise_wilcox_test`, makes it more flexible.

## v0.9.0 - [2025-02-24]

### `Added`
- Created new `preoteomics-class.R`
- Added basic `cli` error-handling in abstract class

## v0.8.3 - [2025-01-14]

### `Fixed`
- Labelling issue of features in `feature_glom`. More robust and equal fast method.

## v0.8.2 - [2025-01-02]

### `Added`
- Finalized docker image for `autoFlow.R`, placed it on docker hub.

## v0.8.1 - [2024-12-17]

### `Changed`
- autoflow now takes optionally pre-computed alpha or beta diversity in table format.

## v0.8.0 - [2024-12-16]

### `Added`
- Added documentation to all functions based on `roxygen2` requirements.
- Created `write_biom` function in sub-class `metataxonomics`

### `Changed`
- Finalised autoFlow.R standalone function, with docker image
- `differential_feature_expression` now uses `sparseMatrix`
- Updated manuals

### `Fixed`
- Error handling in `tools` class
- autoFlow now handles missing data, including all other methods of `tools` class

## v0.7.0 - [2024-11-28]

### `Added`
- Created `hill_taxa.R` that implements `sparseMatrix`.

### `Changed`
- Finalised autoFlow.R standalone function, with docker image
- `differential_feature_expression` now uses `sparseMatrix`

## v0.6.0 - [2024-11-25]

### `Added`
- Created new `diversity` function that uses `sparseMatrix` and uses 5x less memory.

### `Deprecated`
- Removed `vegan::diversity`

## v0.5.0 - [2024-11-14]

### `Changed`
- countData class item now uses `sparseMatrix` instead of `data.table`

## v0.4.1 - [2024-11-8]

### `Added`
- `testthat` R scripts for metataxonomics-class.R

### `Fixed`
- improved `removeZeros` for efficiency.

## v0.4.0 - [2024-10-21]

### `Added`
- Test data
- included documentation, focused on creating an R package

### `Changed`
- alpha_diversity uses `ggpubr` R package
- Improved visualizations in `differential_feature_expression`

## v0.3.1 - [2024-10-13]

### `Added`
- created a `fetch_colors` function to get setNames from `RColorBrewer`.

### `Changed`
- renamed shannon to `alpha_diversity`, more options for shannon metrics

## v0.3.0 - [2024-10-11]

### `Added`
- added dockerfile to use autoFlow

### `Changed`
- Replaced single R scripts for `autoFlow-class.R`
- Applied `foreach` parallel in foldchange computation

### `Deprecated`
- Removed single R scripts for automated analysis

## v0.2.1 - [2024-10-07]

### `Fixed`
- Added private functions within tools class; `original_data`, `tmp_link` and `tmp_restore`. Prevents modifications of tools components.

## v0.2.0 - [2024-09-25]

### `Added`
- Created a new foldchange plot that replaces DFE_plot
- Created `differential_feature_expression` within `tools` class
- Created a volcano_plot

### `Changed`
- Created `tools-class.R`

### `Deprecated`
- Removed deprecated phyloseq functions
- Removed DFE_plot

## v0.1.0 - [2024-09-17]

### `Added`
- Created a tools abstract class based on R6 OOP method, which uses `data.table` for efficiency.
- Created function `removeZeros` within class `tools`.
- Created function `feature_subset` within class `tools`.
- Created function `sample_subset` within class `tools`.
- Created function `feature_glom` within class `tools`.
- Created function `transform` within class `tools`.
- Created function `composition` within class `tools`.
- Created function `shannon` within class `tools`.
- Created function `rankstat` within class `tools`.
- Created function `ordination` within class `tools`.
- Differential Feature Analysis plot; `fold_plot.R`
- Started with `metataxonomics` sub-class for 16S metagenomics data
- Created `print()` similar to Phyloseq
- Created a `reset()` to undo changes
- Created a private `generate_matrix` from a biom file.
- Created an utils folder of new graph functions for `composition_plot`, `ordination_plot`, `paired_fold`, `stats_plot`, `unpaired_fold`
