#!/usr/bin/Rscript

# Load Library -----------------------------------------------------------------
library("OmicFlow")
library("ggplot2")

## Parse args$options from command line
option_list <- list (
  optparse::make_option("--omics",
                        action = "store",
                        type = "character",
                        default = "metagenomics",
                        help="The mode to be selected, default: metagenomics"),
  optparse::make_option(c("-m", "--metadata"),
                        action = "store",
                        type = "character",
                        help="
                        A tab separated file, it should contain a `SAMPLE_ID` column with matching ids in the count table.
                        Optionally, a `SAMPLEPAIR_ID` or `FEATURE_ID` can be provided. Note: spaces are not allowed!"),
  optparse::make_option(c("-b", "--biom"),
                        action = "store",
                        type = "character",
                        help="
                        A Biological Observation Matrix (BIOM) format v2.1.0 compatible to the python module: http://biom-format.org.
                        Supports both BIOM in HDF5 and JSON format."),
  optparse::make_option(c("-t", "--tree"),
                        action = "store",
                        type = "character",
                        help="Phylogenetic tree in newick format, should be supported by `ape::read.tree` (OPTIONAL, default: NULL)",
                        default = NULL),
  optparse::make_option(c("-o", "--outdir"),
                        action = "store",
                        type = "character",
                        help="The directory to write the report, by default the current path is used.",
                        default = getwd()),
  optparse::make_option(c("-f", "--filename"),
                        action = "store",
                        type = "character",
                        help="Name of the HTML report, by default it is named as 'report.html'",
                        default = "report.html"),
  optparse::make_option(c("--i-beta-div"),
                        action = "store",
                        type = "character",
                        help="A custom dissimilarity matrix in TSV/CSV format. Compressed files are also supported."),
  optparse::make_option(c("--i-alpha-div"),
                        action = "store",
                        type = "character",
                        help="A custom alpha diversity created from rarefaction from QIIME2, it should contain columns starting wtih `depth-`"),
  optparse::make_option(c("-c", "--cpus"),
                        action = "store",
                        type = "numeric",
                        help="Number of cores to use, only used when a distance matrix, such as UniFrac or bray-curtis is computed.",
                        default = 4),
  optparse::make_option("--threads",
                        action = "store",
                        type = "numeric",
                        help="Number of threads to be used, mainly for querying of data.tables, default is set at 4",
                        default = 4)
)

parser <- optparse::OptionParser(
  usage = "Usage: Rscript %autoFlow.R [options]",
  option_list = option_list
  )
args <- optparse::parse_args(parser, positional_arguments=TRUE)

## Set threads for data.table
data.table::setDTthreads(threads = args$options$threads)

# main -------------------------------------------------------------------------
# switch statement based on omic selected, create object
if (args$options$omics == "metagenomics") {
  tax <- metagenomics$new(
    metaData = args$options$metadata,
    biomData = args$options$biom,
    treeData = args$options$tree
  )
  tax$feature_subset(Kingdom == "Bacteria")
  tax$normalize()
  
  # Perform automated analysis
  tax$autoFlow(
    beta_div_table = args$options$`i-beta-div`,
    alpha_div_table = args$options$`i-alpha-div`,
    cpus = args$options$cpus,
    normalize = FALSE,
    filename = file.path(args$options$outdir, args$options$filename)
  )
  
} else {
  cli::cli_alert_warning("{args$options$omics} did not match `metagenomics`")
}