#' Wrapper of \link[optparse] for autoFlow.R
#'
#' @description Contains \link[optparse]{make_option} arguments to import required data formats for autoFlow analysis.
#' @return A list of arguments to be read into autoFlow, see \link[optparse]{parse_args}.
#'
#' @export

parse_commandline <- function() {
  option_list <- list (optparse::make_option(c("-m", "--metadata"),
                                             action = "store",
                                             help="tab seperated file"),
                       optparse::make_option(c("-b", "--biom"),
                                             action = "store",
                                             help="biom format file"),
                       optparse::make_option(c("-t", "--tree"),
                                             action = "store",
                                             help="Phylogenetic tree in newick format"),
                       optparse::make_option(c("-c", "--cpus"),
                                             action = "store",
                                             help="Number of cores",
                                             default = 4),
                       optparse::make_option(c("-o", "--outdir"),
                                             action = "store",
                                             help="Output directory",
                                             default = normalizePath(getwd())),
                       optparse::make_option(c("--i-beta-div"),
                                             action = "store",
                                             help="custom beta diversity from qiime2"),
                       optparse::make_option(c("--i-alpha-div"),
                                             action = "store",
                                             help="custom alpha diversity with rarefraction from qiime2")
  )

  parser <- optparse::OptionParser(option_list = option_list)
  arguments <- optparse::parse_args(parser, positional_arguments=TRUE)
  return(arguments$options)
}
