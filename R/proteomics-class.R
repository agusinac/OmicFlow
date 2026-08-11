#' Sub-class proteomics
#'
#' @description This is a sub-class that is compatible to preprocessed data obtained from https://fragpipe.nesvilab.org/. 
#' It inherits all methods from the abstract class \link{omics} and only adapts the \code{initialize} function.
#' It supports pre-existing data structures or paths to text files.
#' When omics data is very large, data loading becomes very expensive. It is therefore recommended to use the [`reset()`](#method-reset) method to reset your changes.
#' Every omics class creates an internal memory efficient back-up of the data, the resetting of changes is an instant process.
#' @seealso \link{omics}
#' @export

proteomics <- R6::R6Class(
  classname = "proteomics",
  cloneable = TRUE,
  inherit = omics,
  active = list(
    #' @field treeData A "phylo" class, see \link[ape]{as.phylo}.
    treeData = function(value) {
      # back-up
      .countData <- private$.countData
      .featureData <- private$.featureData
      .metaData <- private$.metaData
      .treeData <- private$.treeData

      # restore on error
      success <- FALSE
      on.exit({
        if (!success) {
          private$.countData <- .countData
          private$.featureData <- .featureData
          private$.metaData <- .metaData
          private$.treeData <- .treeData
        }
      }, add = TRUE)

      if (missing(value)) {
        success <- TRUE
        private$.treeData
      } else if (inherits(value, "phylo")) {
        private$.treeData <- value
        private$sync()
        self$print()
        success <- TRUE
        invisible(self)
      } else {
        cli::cli_abort("Data input must be {.cls phylo} like {.field treeData}.")
      }
    }
  ),
  public = list(
    #' @description
    #' Initializes the proteomics class object with \code{proteomics$new()}
    #' @param countData A path to an existing file or a dense/sparse \link[Matrix]{Matrix} format.
    #' @param featureData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @param metaData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @param treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree}.
    #' 
    #' @return A new `proteomics` object.
    initialize = function(
      countData = NULL, 
      metaData = NULL, 
      featureData = NULL, 
      treeData = NULL
    ) {
      
      super$initialize(
        countData = countData,
        metaData = metaData,
        featureData = featureData
      )

      # check if `countData` is not empty
      if (is.null(private$.countData))
        cli::cli_abort("{.field countData} cannot be empty.. did you forgot to specify the {.field countData} in {.fun proteomics$new} ?")

      #-------------------#
      ###   treeData    ###
      #-------------------#

      if (!is.null(treeData)) {
        if (is.character(treeData) && length(treeData) == 1 && file.exists(treeData)) {
          private$.treeData <- ape::read.tree(treeData)
          cli::cli_alert_success("{.field treeData} is loaded.")
        } else if (inherits(treeData, "phylo")) {
          private$.treeData <- treeData
          cli::cli_alert_success("{.field treeData} is loaded.")
        } else {
          cli::cli_alert_warning("The provided {.field treeData} could not be loaded. Make sure the tree is supported by {.fun ape::read.tree}")
        }

        # Aligning featureData and countData rows by tree tips
        private$.featureData <- private$.featureData[
          base::order(
            base::match(
              x = private$.featureData[[ private$.feature_id ]], 
              table = private$.treeData$tip.label
            )
          )
        ]
        private$.countData <- private$.countData[private$.featureData[[ private$.feature_id ]], ]
      }
      private$sync()
      self$print()

      # saves data for reset function
      private$original_data = list(
        counts = private$.countData,
        features = private$.featureData,
        metadata = private$.metaData,
        tree = private$.treeData
      )
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
