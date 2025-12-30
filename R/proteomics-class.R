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
    #' @field treeData A "phylo" class, see \link[ape].
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
    #' @param countData A path to an existing file or a dense/sparse \link[Matrix]{Matrix} format.
    #' @param featureData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @param metaData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @param treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree}.
    #' @return A new `proteomics` object.
    initialize = function(countData = NULL, metaData = NULL, featureData = NULL, treeData = NULL) {
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
        private$.featureData <- private$.featureData[order(match(private$.featureData[[ private$.feature_id ]], private$.treeData$tip.label))]
        private$.countData <- private$.countData[private$.featureData[[ private$.feature_id ]], ]
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
    #' Displays parameters of the proteomics class via stdout.
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
      if (length(self$countData) > 0) cat(paste0("## countData:\t[ ", ncol(self$countData), " Samples and ", nrow(self$countData), " Features\t] \n"))
      if (length(self$metaData) > 0) cat(paste0("## metaData:\t[ ", ncol(self$metaData), " Variables and ", nrow(self$metaData), " Samples\t] \n"))
      if (length(self$featureData) > 0) cat(paste0("## featureData:\t[ ", ncol(self$featureData)-1, " Attributes and ", nrow(self$featureData), " Proteins\t] \n"))
      if (length(self$treeData) > 0) cat(paste0("## treeData:\t[ ", length(self$treeData$tip.label), " Tips and ", self$treeData$Nnode, " Nodes\t] \n"))
    }
  ),
  private = list(
    # Private data fields
    #-------------------------#
    .countData = NULL,
    .featureData = NULL,
    .metaData = NULL,
    .treeData = NULL,
    .feature_id = "FEATURE_ID",
    .sample_id = "SAMPLE_ID",
    .samplepair_id = "SAMPLEPAIR_ID",
    original_data = list()
  )
)
