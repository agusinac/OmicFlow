# Load libraries ---------------------------------------------------------------
library("dplyr")
library("tidyr")
library("biomformat")
library("phyloseq")
library("ape")
library("microViz")
library("Biostrings")
library("patchwork")
library("vegan")
library("optparse")
library("rmarkdown")

# Load Functions ---------------------------------------------------------------
# Fetching current path
Rscript <- sub("--file=", "", commandArgs()[4])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))

# Load Data wrangling functions
source(file = paste0(current_path, "99_project-functions.R"))

# Utils
sourceDir(path = paste0(current_path, "utils"))

# Parse command line
data_00 <- parse_commandline()

# Creates required directories
create.dir("results", showWarnings = FALSE)
create.dir("RDS", showWarnings = FALSE)

# Run all scripts --------------------------------------------------------------
source(file = paste0(current_path, "01_load.R"))
source(file = paste0(current_path, "02_clean.R"))
source(file = paste0(current_path, "03_data_visualization.R"))
source(file = paste0(current_path, "04_Model-UniFrac-PCoA.R"))
source(file = paste0(current_path, "05_Model-RDA.R"))
rmarkdown::render(input = paste0(current_path, "../documents/Report.Rmd"),
                  output_file = "trial_report.html",
                  output_dir = "results")
