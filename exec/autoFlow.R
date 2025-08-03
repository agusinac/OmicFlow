#!/usr/bin/Rscript

# Load Library -----------------------------------------------------------------
library("OmicFLow")

# Parse command line
options <- parse_commandline()

# main -------------------------------------------------------------------------
# switch statement based on omic selected, create object
omics <- switch(
    options$omics,
    "metagenomics" = metagenomics$new(
        metaData = options$metadata,
        biomData = options$biom,
        treeData = ifelse(!is.null(options$tree), options$tree, NA)
        )
)

# Perform automated analysis
omics$autoFlow(
    beta_div_table = options$`i-beta-div`,
    alpha_div_table = options$`i-alpha-div`,
    cpus = options$cpus,
    filename = paste0(options$outdir, "/", options$filename)
)