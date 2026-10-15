# OmicFlow

## Overview

OmicFlow is a generalised data structure for fast and efficient loading
of various sparse `omics` data. The main `omics` class extends to
`metagenomics` data in text or
[BIOM](https://biom-format.org/documentation/format_versions/biom-2.0.html)
and to `proteomics`. It also supports non-sparse data, but it’s
performance peaks in sparsity.

## Installation

The latest stable version can be installed from CRAN.

``` r

install.packages('OmicFlow', dependencies = TRUE)
```

The development version is available on GitHub.

``` r

install.packages('pak') # if not yet installed
pak::pkg_install('agusinac/OmicFlow@dev')
```

## Usage

Initialize the `metagenomics` or any `omics` object from a file-path or
pre-loaded object.

``` r
library("OmicFlow")

## Load example data
metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")
biom_file <- system.file("extdata", "biomv2.biom", package = "OmicFlow")

## Initialise metagenomics from biom_file
taxa <- metagenomics$new(
    metaData = metadata_file,
    biomData = biom_file,
    treeData = tree_file            # optional
)

## Initialise omics class
taxa2 <- omics$new(
    metaData = metadata_file,
    ## counts_file is allowed to have non-sample ids
    countData = counts_file,
    featureData = features_file,    # optional
    treeData = tree_file            # optional
)

taxa$feature_subset(Kingdom == "Bacteria")
taxa$scale(method = "tss")

# Access variables directly
taxa$metaData
taxa$countData
taxa$featureData
taxa$treeData

# Change variables & enjoy the automatic sync
taxa$featureData <- taxa$featureData[1:100, ]

# Inspect what functions variables are available to the class
taxa$
```

If you are new to `OmicFlow`, the best place to start is the
[Introduction to
OmicFlow](https://agusinac.github.io/OmicFlow/articles/getting-started.html).

## Docker

> \[!NOTE\] `agusinac/autoflow:latest` will be superseded by
> `agusinac/omicflow:latest` from v1.7.0

Example: **Outputs a `report.html` file in current work directory**

``` bash
docker pull agusinac/autoflow:latest

docker run -it --rm -v \
    "$(pwd)":/data \             # Mount the data in a temporary directory
    -w /data \                   # set working directory
    -u $(id -u):$(id -g) \       # non-root user
    agusinac/autoflow:latest \
    autoflow \                   # autoflow R script
    -b /data/biom_with_taxonomy_hdf5.biom \
    -m /data/metadata.tsv
```

## Support

If you are having issues, please [create a
ticket](https://github.com/agusinac/OmicFlow/issues)

## Acknowledgements

- Sticker logo created by [studio
  floorescent](https://www.floorescent.nl/) \| [Floor
  Baas](https://nl.linkedin.com/in/floorbaas)
