#' @title Sub-class proteomics
#' @docType class
#' @section Introduction:
#' OmicFlow is constructed around the main superclass \link{omics} that in turn is 
#' inherited by \link{metagenomics} and \link{proteomics}, that offer additional 
#' fields or functions. The superclass \link{omics} contains both public and private 
#' methods, only the public componenents are documented, whereas the private 
#' methods can be accessed via `proteomics$private_methods`, which is not supported.
#' 
#' The \link{omics} class components use \link{Matrix} and \link{data.table} in the 
#' background for fast loading and data wrangling, these data structures return 
#' the data by reference, so does the R6 class, therefore classes in OmicFlow do 
#' not create copies of the object unless called via \code{$copy()}.
#' 
#' All classes needs to be first initialised via the \code{$new()} method and 
#' require the `metaData` and `countData` field components. Alternatively, you 
#' can also initialise with only the `metaData` and later add the other fields 
#' via the active binding, but this will not create a back-up and is not advised!
#' For more hands-on examples see [Getting Started with OmicFlow](https://agusinac.github.io/OmicFlow/articles/getting-started.html).
#' 
#' \subsection{Additional arguments in \code{proteomics$new()}}{
#' The \code{proteomics} class also offers the loading of a phylogenetic tree as
#' a `phylo` class (see \link[ape]{as.phylo}) that can be supplied to `treeData`. 
#' When a `treeData` is supplied than every class field is arranged according to 
#' the tree tip labels.
#' }
#' 
#' @inheritSection omics Input requirements
#' @inheritSection omics Metadata validation
#' @export

proteomics <- R6::R6Class(
  classname = "proteomics",
  cloneable = FALSE,
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
        cli::cli_abort("Input must be {.cls phylo} like {.field treeData}.")
      }
    }
  ),
  public = list(
    #' @param treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree} (default: \code{NULL}).
    #' @examples
    #' library("OmicFlow")
    #'
    #' ## Method 1: load from filepath
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' obj <- proteomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file
    #' )
    #' 
    #' ## Method 2: Load from data.frame or matrix
    #' n_cols <- 5
    #' n_rows <- 100
    #' n_vals <- n_cols * n_rows
    #' metadata <- data.frame("SAMPLE_ID" = paste0("Sample_", 1:n_cols))
    #' features <- data.frame("FEATURE_ID" = paste0("protein_", 1:n_rows))
    #' counts <- Matrix::Matrix(
    #'  1:n_vals, nrow = n_rows, ncol = n_cols, 
    #'  dimnames = list(features$FEATURE_ID, metadata$SAMPLE_ID)
    #' )
    #' 
    #' obj <- proteomics$new(
    #'  metaData = metadata,
    #'  featureData = features, # optional
    #'  countData = counts
    #' )
    #'
    #' ## Method 3: you have features and counts in a single data.frame
    #' metadata <- data.frame(SAMPLE_ID = c("S1", "S2", "S3"))
    #' counts <- data.frame(S1 = c(2,3,0), S2 = c(2,"",1), S3 = c(2,1,NA), proteins = c("ZEB1", "ZEB2", "MAPK"))
    #' 
    #' obj <- proteomics$new(
    #'  metaData = metadata,
    #'  countData = counts
    #' )
    #' @return A new \link{proteomics} object.
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
