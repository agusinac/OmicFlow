#' Sub-class proteomics
#'
#' @description This is a sub-class that is compatible to preprocessed data obtained from https://fragpipe.nesvilab.org/. 
#' It inherits all methods from the abstract class \link{omics} and only adapts the \code{initialize} function.
#' It supports pre-existing data structures or paths to text files.
#' When omics data is very large, data loading becomes very expensive. It is therefore recommended to use the [`reset()`](#method-reset) method to reset your changes.
#' Every omics class creates an internal memory efficient back-up of the data, the resetting of changes is an instant process.
#' @seealso \link{omics}
#' @import R6
#' @importFrom ape read.tree
#' @export

proteomics <- R6::R6Class(
  classname = "proteomics",
  cloneable = TRUE,
  inherit = omics,
  active = list(
    treeData = function(value) {
      # Restores omics class components
      private$tmp_link(
        .countData = private$.countData,
        .featureData = private$.featureData,
        .metaData = private$.metaData,
        .treeData = private$.treeData
      )
      
      # Returns on failure
      success <- FALSE
      on.exit({
        if (!success) {
          private$tmp_restore()
        }
      }, add = TRUE)

      if (missing(value)) {
        private$.treeData
      } else if (inherits(value, "phylo")){
        private$.treeData <- value
        private$sync()
        success <- TRUE
        invisible(self)
      } else stop("Data input requires to be of the same class as `treeData`")
    }
  ),
  public = list(
    #' @description
    #' Initializes the proteomics class object with \code{proteomics$new()}
    #' @param countData A path to an existing file, data.table, data.frame, matrix or sparseMatrix with zero values.
    #' @param featureData A path to an existing file, data.table or data.frame.
    #' @param metaData A path to an existing file, data.table or data.frame.
    #' @param treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree}.
    #' @return A new `proteomics` object.
    initialize = function(countData = NA, metaData = NA, featureData = NA, treeData = NA) {
      super$initialize(countData = countData,
                       metaData = metaData,
                       featureData = featureData)

      #-------------------#
      ###   treeData    ###
      #-------------------#

      if (!is.null(treeData)) {
        if (is.character(treeData) && length(treeData) == 1 && file.exists(treeData)) {
          private$.treeData <- ape::read.tree(treeData)
          cli::cli_alert_success("treeData is loaded.")
        } else if (inherits(treeData, "phylo")) {
          private$.treeData <- treeData
          cli::cli_alert_success("treeData is loaded.")
        } else {
          cli::cli_alert_warning("The provided TreeData could not be loaded. Make sure the tree is supported by `ape::read.tree`")
        }

        # Aligning featureData and countData rows by tree tips
        private$.featureData <- private$.featureData[order(match(private$.featureData[[ self$feature_id ]], private$.treeData$tip.label))]
        private$.countData <- private$.countData[private$.featureData[[ self$feature_id ]], ]
      }

      self$print()

      # saves data for reset function
      private$original_data = list(
        counts = private$.countData,
        features = private$.featureData,
        metadata = private$.metaData,
        tree = private$.treeData
      )
    },
    #' @description
    #' Displays parameters of the proteomics object via stdout.
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #' tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")
    #'
    #' prot <- proteomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #'  treeData = tree_file
    #' )
    #'
    #' # method 1 to call print function
    #' prot
    #'
    #' # method 2 to call print function
    #' prot$print()
    #'
    #' @return object in place
    print = function() {
      cat("## proteomics-class object \n")
      if (length(private$.countData) > 0) cat(paste0("## countData:\t[ ", ncol(private$.countData), " Samples and ", nrow(private$.countData), " Features\t] \n"))
      if (length(private$.metaData) > 0) cat(paste0("## metaData:\t[ ", ncol(private$.metaData), " Variables and ", nrow(private$.metaData), " Samples\t] \n"))
      if (length(private$.featureData) > 0) cat(paste0("## featureData:\t[ ", ncol(private$.featureData)-1, " Attributes and ", nrow(private$.featureData), " Proteins\t] \n"))
      if (length(private$.treeData) > 0) cat(paste0("## treeData:\t[ ", length(private$.treeData$tip.label), " Tips and ", private$.treeData$Nnode, " Nodes\t] \n"))
    },
        #' @description
    #' Upon creation of a new `proteomics` object a small backup of the original data is created.
    #' Since modification of the object is done by reference and duplicates are not made, it is possible to `reset` changes to the class.
    #' The methods from the abstract class \link{omics} also contains a private method to prevent any changes to the original object when using methods such as \code{ordination} \code{alpha_diversity} or \code{$DFE}.  
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #' tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")
    #'
    #' prot <- proteomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #'  treeData = tree_file
    #' )
    #' 
    #' # Performs modifications
    #' prot$transform(log2)
    #'
    #' # resets
    #' prot$reset()
    #'
    #' @return object in place
    reset = function() {
      private$.countData = private$original_data$counts
      private$.featureData = private$original_data$features
      private$.metaData = private$original_data$metadata
      private$.treeData = private$original_data$tree
      invisible(self)
    }
  ),
  private = list(
    # Private data fields
    #-------------------------#
    .countData = NULL,
    .featureData = NULL,
    .metaData = NULL,
    .treeData = NULL,
    original_data = list()
  )
)
