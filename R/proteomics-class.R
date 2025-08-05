#' Sub-class proteomics
#'
#' @description This is a sub-class for proteomics data, called proteomics.
#' It inherits all methods from the abstract class \link{omics} and only adapts the \code{initialize} function.
#'
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

      if (!is.na(treeData)) self$treeData <- ape::read.tree(treeData)

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
    #' Displays parameters of the `proteomics` object via stdout.
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
    #' The methods from the abstract class `omics` also contain a private method to prevent any changes to the original object. Such cases are ordination, alpha_diversity, differential_feature_expression.
    reset = function() {
      self$countData = private$original_data$counts
      self$featureData = private$original_data$features
      self$metaData = private$original_data$metadata
      self$treeData = private$original_data$tree
      invisible(self)
    },
    #' @description
    #' Removes empty (zero) values by row, column and tips.
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
