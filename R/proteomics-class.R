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
  cloneable = FALSE,
  inherit = omics,
  public = list(
    #' @field countData A path to an existing file, data.table or data.frame.
    countData = NULL,

    #' @field metaData A path to an existing file, data.table or data.frame.
    metaData = NULL,

    #' @field featureData A path to an existing file, data.table or data.frame.
    featureData = NULL,

    #' @field treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree}.
    treeData = NULL,

    #' @description
    #' Initializes the proteomics class object with \code{proteomics$new()}
    #' @param countData countData A path to an existing file or sparseMatrix.
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
          self$treeData <- ape::read.tree(treeData)
          cli::cli_alert_success("treeData is loaded.")
        } else if (inherits(treeData, "phylo")) {
          self$treeData <- treeData
          cli::cli_alert_success("treeData is loaded.")
        } else {
          cli::cli_alert_warning("The provided TreeData could not be loaded. Make sure the tree is supported by `ape::read.tree`")
        }
      }

      self$print()

      # saves data for reset function
      private$original_data = list(
        counts = self$countData,
        features = self$featureData,
        metadata = self$metaData,
        tree = self$treeData
      )
    },
    #' @description
    #' Displays parameters of the proteomics object via stdout.
    #' @examples
    #' prot_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' prot <- readRDS(prot_path)
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
    },
        #' @description
    #' Upon creation of a new `proteomics` object a small backup of the original data is created.
    #' Since modification of the object is done by reference and duplicates are not made, it is possible to `reset` changes to the class.
    #' The methods from the abstract class \link{omics} also contains a private method to prevent any changes to the original object when using methods such as \code{ordination} \code{alpha_diversity} or \code{$DFE}.  
    #' @examples
    #'  
    #' prot_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' prot <- readRDS(prot_path)
    #' 
    #' # Performs modifications
    #' prot$transform(log2)
    #'
    #' # resets
    #' prot$reset()
    #'
    #' @return object in place
    reset = function() {
      self$countData = private$original_data$counts
      self$featureData = private$original_data$features
      self$metaData = private$original_data$metadata
      self$treeData = private$original_data$tree
      invisible(self)
    },
    #' @description
    #' Removes empty (zero) values by row, column and tips from the `countData` and `treeData`.
    #' This method is performed automatically during subsetting of the object.
    #' @importFrom ape keep.tip
    #' @examples
    #' prot_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' prot <- readRDS(prot_path)
    #' 
    #' # Sample subset induces empty features
    #' prot$sample_subset(treatment == "tumor")
    #'
    #' # Remove empty features from countData and treeData
    #' prot$removeZeros()
    #' 
    #' @return object in place
    removeZeros = function() {
      super$removeZeros()
      if (!is.null(self$treeData)) self$treeData <- ape::keep.tip(self$treeData, self$featureData$FEATURE_ID)
      invisible(self)
    }
  ),
  private = list(
    original_data = list()
  )
)
