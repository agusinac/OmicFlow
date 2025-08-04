#!/usr/bin/Rscript

# Load Library -----------------------------------------------------------------
library("OmicFLow")

# Parse command line
options <- parse_commandline()

# main -------------------------------------------------------------------------
# switch statement based on omic selected, create object
if (options$omics == "metagenomics") {
  tax <- metagenomics$new(
    metaData = options$metadata,
    biomData = options$biom,
    treeData = ifelse(!is.null(options$tree), options$tree, NULL)
  )
  tax$feature_subset(Kingdom == "Bacteria")
  tax$normalize()
  
  # Perform automated analysis
  tax$autoFlow(
    beta_div_table = options$`i-beta-div`,
    alpha_div_table = options$`i-alpha-div`,
    cpus = options$cpus,
    normalize = FALSE,
    filename = paste0(options$outdir, "/", options$filename)
  )
  
} else {
  cli::cli_alert_warning("{options$omics} did not match `metagenomics`")
}