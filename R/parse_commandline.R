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
                       optparse::make_option(c("-o", "--outdir"),
                                   action = "store",
                                   help="Output directory")
  )
  
  parser <- optparse::OptionParser(option_list = option_list)
  arguments <- optparse::parse_args(parser, positional_arguments=TRUE)
  return(arguments$options)
}