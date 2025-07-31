# Load Library -----------------------------------------------------------------
library("Matrix")
library("ggplot2")
library("ggtree")
library("patchwork")

# Fetch current Rscript path
Rscript <- sub("--file=", "", grep("^--file", commandArgs(), value = TRUE)[1])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))
# Since OmicFlow is not a package yet, current work-around:
devtools::load_all(current_path)

# Parse command line
options <- parse_commandline()

outfile_path <- options$outdir

# main -------------------------------------------------------------------------
# switch statement based on omic selected, create object
omics <- metataxonomics$new(metaData = options$metadata,
                            biomData = options$biom,
                            treeData = ifelse(!is.null(options$tree), options$tree, NA))

# Set parameters
feature_ranks = c("Phylum", "Family", "Genus")
distance_metrics = c("unifrac")
feature_filter = c("uncultured")
rankstat_labels <- sub("RANKSTAT_", "", colnames(omics$metaData)[grepl("RANKSTAT_", colnames(omics$metaData))])

# Fetch additional files from preprocessing pipeline ---------------------------
# file_location <- normalizePath(getwd())
# read_stats_path <- paste0(file_location, "/read_stats.txt")
# sankeyplot <- paste0(file_location, "/sankeyplots/biotaviz_sankey_prepfile-AverageAllSamples.png")
#
# if (file.exists(read_stats_path)) {
#   preprocess_stats <- readr::read_delim(read_stats_path)
# } else {
#   preprocess_stats <- NULL
# }
#
# if (file.exists(sankeyplot)) {
#   sankeyplot <- sankeyplot
# } else {
#   sankeyplot <- NULL
# }

# Initiate OmicFlow
plots <- omics$autoFlow(feature_ranks = feature_ranks,
                       feature_filter = feature_filter,
                       distance_metrics = distance_metrics,
                       dist_matrix = options$`i-beta-div`,
                       alpha_div_table = options$`i-alpha-div`,
                       cpus = options$cpus)

# Create report
rmarkdown::render(input = paste0(current_path, "report/report.Rmd"),
                  output_file = paste0(outfile_path, "/report.html"))
