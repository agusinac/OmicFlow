# Load libraries ---------------------------------------------------------------
library("dplyr")
library("ggplot2")
library("foreach")
library("patchwork")
library("optparse")

# Load Functions ---------------------------------------------------------------
sourceDir <- function(path, trace = TRUE, ...) {
  op <- base::options(); on.exit(base::options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}
# Fetch current Rscript path
Rscript <- sub("--file=", "", commandArgs()[4])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))

# Data wrangling functions
sourceDir(path = paste0(current_path, "R"))

# Parse command line
data_00 <- parse_commandline()
outfile_path <- normalizePath(data_00$outDir)

# Taxa to be visualized:
taxa_names <- c("Phylum", "Family", "Genus")

# main -------------------------------------------------------------------------
# Create new object class

# Run autoFlow

# pass autoFlow plots list to report.rmd (should automatically detect number of lists and their formats)
rmarkdown::render(input = paste0(getwd(), "/automated-omics-analysis/documents/report.Rmd"),
                  output_file = paste0("report.html"))
