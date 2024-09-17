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
                       optparse::make_option(c("-r", "--refseq"),
                                   action = "store",
                                   help="Reference sequence in fasta format"),
                       optparse::make_option(c("-o", "--outdir"),
                                   action = "store",
                                   help="Output directory")
  )
  
  parser <- optparse::OptionParser(option_list = option_list)
  arguments <- optparse::parse_args(parser, positional_arguments=TRUE)
  
  opt <- arguments$options
  print(opt)
  if (base::grepl(".tsv$", opt$metadata)) {
    metadata <- phyloseq::import_qiime_sample_data(mapfilename = opt$metadata)
  } else stop("Please provide a tab-separated metadata file!")
  
  if (base::grepl(".biom$", opt$biom)) {
    biom_data <- phyloseq::import_biom(BIOMfilename = opt$biom)
  } else stop("Please provide a valid biom file")

  if (!is.null(opt$tree)) {
    if (base::grepl(".newick$", opt$tree)) {
      tree <- phyloseq::read_tree(treefile = opt$tree)
    }
  } else {
    tree <- NULL
  }

  if (!is.null(opt$refseq)) {
    if (base::grepl(".fasta$", opt$refseq)) {
      refseq <- Biostrings::readDNAStringSet(filepath = opt$refseq)
    }
  } else {
    refseq <- NULL
  }

  result <- list(
    metaData = metadata,
    biomData = biom_data,
    treeData = tree,
    refseqData = refseq,
    outDir = opt$outdir
  )

  return(result)
}