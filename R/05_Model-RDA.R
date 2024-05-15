# Parsing command line --------------------------------------------------------#
script_parser <- function() {
  option_list <- list (optparse::make_option(c("--inomics"), 
                                             action = "store",
                                             default = NULL,
                                             help="sample counts .xlsx file"),
                       optparse::make_option(c("--osheet"),
                                             action = "store",
                                             default = NULL,
                                             help="excel sheet for omics"),
                       optparse::make_option(c("--inpheno"),
                                             action = "store",
                                             default = NULL,
                                             help="taxonomical .xlsx file"),
                       optparse::make_option(c("--psheet"),
                                             action = "store",
                                             default = NULL,
                                             help="excel sheet for taxonomical file"),
                       optparse::make_option(c("--contrast"),
                                             action = "store",
                                             default = NULL,
                                             help="Metadata column to be compared against sample counts"),
                       optparse::make_option(c("--log"),
                                             action = "store",
                                             default = 1,
                                             help="Log transformation base"),
                       optparse::make_option(c("--output"),
                                             action = "store",
                                             default = ".",
                                             help="Output directory"),
                       optparse::make_option(c("--pairwise"),
                                             action = "store",
                                             default = FALSE,
                                             help = "Boolean operator, TRUE or FALSE, default = FALSE")
  )
  
  parser <- optparse::OptionParser(option_list = option_list)
  arguments <- optparse::parse_args(parser, positional_arguments=TRUE)
  
  opt <- arguments$options
  print(opt)
  # Import omics df
  if (!is.null(opt$inomics)) {
    if (base::grepl(".xlsx$", opt$inomics) && opt$osheet != "") {
      # imports as a Tibble, but these don't allow for row names [ https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html ]
      omics <- readxl::read_excel(opt$inomics, sheet = opt$osheet)
      omics <- as.data.frame(omics)
      base::rownames(omics) <- omics[,1]
      omics[,1] <- NULL
    }
  }
  
  # Import metadata df
  if (!is.null(opt$inpheno)) {
    if (base::grepl(".xlsx$", opt$inpheno) && opt$psheet != "") {
      # imports as a Tibble, but these don't allow for row names [ https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html ]
      meta <- readxl::read_excel(opt$inpheno, sheet = opt$psheet)
      meta <- as.data.frame(meta)
      rownames(meta) <- meta[,1]
      meta[,1] <- NULL
    }
  }
  
  if (!is.null(opt$contrast)) {
    contrast <- opt$contrast
  }
  
  if (!is.null(opt$inomics) && !is.null(opt$inpheno)) {
    result <- list(
      otu_tab = omics,
      meta_tab = meta,
      contrast = contrast,
      log = opt$log,
      outdir = opt$output,
      pairwise = opt$pairwise
    )
  } else {
    result <- NULL
  }
  return(result)
}

logn <- function(otu_tab, scalar=1) {
  # log-transform, center
  # Y' = log ( A * Y + 1 ) ; where A is the 'strength' of the log transformation : 1, 10, 100, 1000, etc., default = 1
  otu_tab.log <- ( scalar * otu_tab ) + 1 
  otu_tab.log <- log( otu_tab.log ) 
  otu_tab.sc <- scale(otu_tab.log, center = TRUE, scale = FALSE)
  return(otu_tab.sc)
}

# Main code -------------------------------------------------------------------#
# Checks if arguments are passed from 05_main.R
if (grepl("05_Model-RDA.R", commandArgs()[4])) {
  data_05 <- script_parser()
  otu_tab <- data_05$otu_tab
  meta_tab <- data_05$meta_tab
  contrast <- data_05$contrast
  
  # Loads required functions/libraries
  source("utils/pairwise_triplot.R")
  
  # Models RDA
  otu_tab.log <- logn(otu_tab, scalar = data_05$log)
  mod.rda <- vegan::rda(otu_tab.log ~ get(contrast, meta_tab) + Condition(NULL), 
                        data = meta_tab, 
                        scale = FALSE, 
                        na.action = na.fail, 
                        subset = NULL)
  
  if (data_05$pairwise) {
    mod.pca <- vegan::rda(otu_tab.log,
                          scale = FALSE)
    
    pairwise_triplot(model = mod.pca,
                     target_col = contrast,
                     metadata = meta_tab,
                     pairwise = data_05$pairwise,
                     outdir = data_05$outdir)
  }
  
  # Saves plot
  pdf(file=paste0(data_05$outdir, "/", "RDA_triplot.pdf"))
  pairwise_triplot(model = mod.rda, 
                   target_col = contrast, 
                   metadata = meta_tab,
                   choice_dim = c("RDA1", "PC1"))
  pairwise_triplot(model = vegan::rda(otu_tab.log, scale=FALSE), 
                   target_col = contrast, 
                   metadata = meta_tab,
                   choice_dim = c("PC1", "PC2"))
  dev.off()
}
# Checks if arguments are passed from 00_main.R
if (grepl("00_main.R", commandArgs()[4])) {
  radboud_ps_rel_bac_norm <- readRDS("data/RDS/02_ps_rel_bac_norm.rds")
  otu_tab <- get_otu(ps = radboud_ps_rel_bac_norm)
  meta_tab <- get_meta(ps = radboud_ps_rel_bac_norm)
  
  # Natural log transformation + scaling
  otu_tab.log <- logn(otu_tab, scalar = 10)
  # Model selection: Outputs significant columns
  mod.res <- rda_select(otu_tab.log, meta_table = meta_tab)
  
  # Fetch significant columns
  mod.cols <- mod.res["model.col"]
  
  # PCA model
  mod.pca <- vegan::rda(otu_tab.log,
                        scale = FALSE)
  
  # Create empty plot list
  plot_list <- list()
  n = 0
  # Creates triplot of different columns
  for (i in mod.cols) {
    # RDA model
    mod.rda <- vegan::rda(otu_tab.log ~ get(i, meta_tab) + Condition(NULL),
                          data = meta_tab,
                          scale = FALSE,
                          na.action = na.fail,
                          subset = NULL)
    n = n + 1
    triplot_rda <- function(){pairwise_triplot(model = mod.rda,
                                               target_col = i,
                                               metadata = meta_tab,
                                               pairwise = FALSE,
                                               choice_dim = c("RDA1", "PC1"))
    }
    plot_list[[n]] <- triplot_rda
    n = n + 1
    triplot_pca <- function(){pairwise_triplot(model = mod.pca,
                                               target_col = i,
                                               metadata = meta_tab,
                                               pairwise = FALSE,
                                               choice_dim = c("PC1", "PC2"))
    }
    plot_list[[n]] <- triplot_pca
  }
  
}
