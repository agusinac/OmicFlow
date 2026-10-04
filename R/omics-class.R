#' @title Superclass omics
#' @docType class
#' @section Introduction:
#' OmicFlow is constructed around the main superclass `omics` that in turn is 
#' inherited by \link{metagenomics} and \link{proteomics}, that offer additional 
#' fields or functions. The superclass `omics` contains both public and private 
#' methods, only the public componenents are documented, whereas the private 
#' methods can be accessed via `omics$private_methods`, which is not supported.
#' 
#' The `omics` class components use \link{Matrix} and \link{data.table} in the 
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
#' @section Input requirements:
#' Every class in OmicFlow requires at least a `metaData` and `countData` field to
#' be specified. The `metaData` must contain a `SAMPLE_ID` and optionally a `SAMPLEPAIR_ID`, 
#' see section "Metadata validation" below.
#' 
#' The `countData` can be both a matrix or data frame. If the original contains 
#' rownames, then these will become the `FEATURE_ID` in `featureData` (automatically created), 
#' otherwise a random feature identifier is assigned. If the `countData` is a filetype
#'  or a data frame with both matching sample identifiers in `SAMPLE_ID` and other 
#' additional columns, then these will be automatically split, columns matching the 
#' `SAMPLE_ID` will be placed in the `countData`, any other columns will be placed 
#' in front of the `FEATURE_ID` column in `featureData`.
#' 
#' @section Metadata validation:
#' Every class in OmicFlow begins by validating the sample metadata. 
#' This validation checks that sample identifiers match the sample IDs
#' present in the abundance tables.
#'
#' The validation rules are defined in JSON format in the omics class. 
#' Metadata can be supplied as either:
#' \itemize{
#'   \item A CSV or TSV file.
#'   \item A `data.table`.
#' }
#'
#' In both cases, metadata must include a header row and follow a
#' one-row-per-sample structure. For file inputs, the header must be 
#' on the first line. Columns not described below are permitted and 
#' are ignored during metadata validation.
#'
#' \subsection{Minimum requirement}{
#' The metadata must contain a `SAMPLE_ID` column.
#'
#' \itemize{
#'   \item Each row must have a unique, non-empty sample identifier.
#'   \item Sample identifiers cannot contain spaces; use underscores (`_`) or
#'   dashes (`-`) instead.
#' }
#'
#' Example metadata:
#'
#' \preformatted{
#' SAMPLE_ID  SAMPLEPAIR_ID  CONTRAST_Treatment  VARIABLE_Age
#' S1         P1             Drug                42
#' S2         P1             Placebo             36
#' S3         P2             Drug                51
#' }
#' }
#'
#' \subsection{Column types and naming rules}{
#'
#' \strong{Required column}
#'
#' \tabular{lll}{
#' Column \tab Type \tab Rules \cr
#' `SAMPLE_ID` \tab string \tab Unique, non-empty, no spaces; one identifier
#' per sample row
#' }
#'
#' \strong{Optional standard columns}
#'
#' \tabular{lll}{
#' Column \tab Type \tab Rules \cr
#' `SAMPLEPAIR_ID` \tab string \tab Optional; no spaces. Used for samples that
#' are paired or originate from the same individual or subject.
#' }
#' Some functions require a `SAMEPLEPAIR_ID` in the `metaData` when \code{paired = TRUE} is used.
#' The user can add or remove the `SAMPLEPAIR_ID` at anytime to the `metaData` via active binding.
#'
#' \strong{Pattern-based columns}
#'
#' Additional variables can be defined with special column-name prefixes:
#'
#' \itemize{
#'   \item `CONTRAST_...`: Grouping or category labels for differential
#'   comparisons. For example, `CONTRAST_Treatment` may contain `Drug` and
#'   `Placebo`.
#'   \item `VARIABLE_...`: Numeric or string variables for statistical
#'   analysis. For example, `VARIABLE_Age` may contain values such as `42` and
#'   `51`.
#' }
#' }
#' Pattern-based columns are only used in [`autoFlow()`](#method-autoFlow). Currently, only columns
#' with the `CONTRAST_` prefix are supported.
#' 
#' @param countData A path to an existing file, \link[Matrix]{Matrix}, \link[data.table]{data.table} or \link[base]{data.frame} (default: \code{NULL}).
#' @param featureData A path to an existing file, \link[data.table]{data.table} or \link[base]{data.frame} (default: \code{NULL}).
#' @param metaData A path to an existing file, \link[data.table]{data.table} or \link[base]{data.frame} (default: \code{NULL}).
#' @export

omics <- R6::R6Class(
  classname = "omics",
  cloneable = TRUE,
  active = list(
    #' @field metaData A \link[data.table]{data.table} with `SAMPLE_ID` column.
    metaData = function(value) {
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
        private$.metaData
      } else if (inherits(value, "data.table")) {
        private$.metaData <- value
        private$sync()
        success <- TRUE
        self$print()
        invisible(self)
      } else {
        cli::cli_abort("Input must be {.cls data.table} like {.field metaData}.")
      }
    },
    #' @field featureData A \link[data.table]{data.table} with `FEATURE_ID` column.
    featureData = function(value) {
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
        private$.featureData
      } else if (inherits(value, "data.table")) {
        private$.featureData <- value
        private$sync()
        success <- TRUE
        self$print()
        invisible(self)
      } else {
        cli::cli_abort("Input must be {.cls data.table} like {.field featureData}.")
      }
    },
    #' @field countData A \link[Matrix]{sparseMatrix} with rownames matching `FEATURE_ID` of `featureData` and colnames matching `SAMPLE_ID` of `countData`.
    countData = function(value) {
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
        private$.countData
      } else if (inherits(value, "Matrix")) {
        private$.countData <- value
        private$sync()
        success <- TRUE
        self$print()
        invisible(self)
      } else {
        cli::cli_abort("Input must be {.cls sparseMatrix} like {.field countData}.")
      }
    }
  ),
  public = list(
    #' @examples
    #' library("OmicFlow")
    #'
    #' ## Method 1: load from filepath
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' obj <- omics$new(
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
    #' obj <- omics$new(
    #'  metaData = metadata,
    #'  featureData = features, # optional
    #'  countData = counts
    #' )
    #'
    #' ## Method 3: you have features and counts in a single data.frame
    #' metadata <- data.frame(SAMPLE_ID = c("S1", "S2", "S3"))
    #' counts <- data.frame(S1 = c(2,3,0), S2 = c(2,"",1), S3 = c(2,1,NA), proteins = c("ZEB1", "ZEB2", "MAPK"))
    #' 
    #' obj <- omics$new(
    #'  metaData = metadata,
    #'  countData = counts
    #' )
    #' @return A new \link{omics} object.
    initialize = function(
      metaData = NULL,
      countData = NULL, 
      featureData = NULL
    ) {
      #-------------------#
      ###   metaData    ###
      #-------------------#
      if (!is.null(metaData)) {
        duplicated_sample_ids <- FALSE

        ## check if message is returned
        private$.metaData <- private$check_table(metaData)
        if (!data.table::is.data.table(private$.metaData)) {
          cli::cli_abort(c(
            "Error in {.field metaData}:",
            "x" = cli::format_inline("{private$.metaData}")
          ))
        }
        private$validate()

        if (private$.valid_schema) {
          cli::cli_alert_success("{.field metaData} template passed the JSON validation.")

          #--------------------------------------------------------------------#
          ## Checking for duplicated sample identifiers
          #--------------------------------------------------------------------#

          cli::cli_alert_info("Checking for duplicated identifiers ..")
          duplicated_sample_idx <- base::duplicated(private$.metaData, by = private$.sample_id)
          duplicated_sample_ids <- any(duplicated_sample_idx)
          if (duplicated_sample_ids) {
            duplicated_sample_names <- data.table::unique(private$.metaData[[private$.sample_id]][duplicated_sample_idx])
            cli::cli_abort(
              "Found duplicated: {.val {duplicated_sample_names}}\
              \n Make sure {.arg SAMPLE_ID} column contains {.strong unique} identifiers!"
            )
          }

        } else {
          errors <- attr(private$.valid_schema, "errors")
          cli::cli_abort(
            "JSON validation failed: \n{ paste(errors$message, collapse = '\n')}"
            )
        }

      } else {
        cli::cli_abort(
          "{.field metaData} cannot be empty, please provide a {.cls data.frame}, {.cls data.table} or {.val filepath}"
        )
      }

      #-------------------#
      ###  featureData  ###
      #-------------------#
      if (!is.null(featureData)) {
        duplicated_feature_ids <- FALSE

        private$.featureData <- private$check_table(featureData)
        if (!data.table::is.data.table(private$.featureData)) {
          cli::cli_abort(c(
            "Error in {.field featureData}:",
            "x" = cli::format_inline("{private$.featureData}")
          ))
        }
        
        if (column_exists(private$.feature_id, private$.featureData)) {
          duplicated_feature_idx <- base::duplicated(private$.featureData, by = private$.feature_id)
          duplicated_feature_ids <- any(duplicated_feature_idx)

          if (duplicated_feature_ids) {
            duplicated_feature_names <- data.table::unique(private$.featureData[[private$.feature_id]][duplicated_feature_idx])
            cli::cli_abort(
              "Found duplicated: {.val {duplicated_feature_names}} \
              \n Make sure {.arg FEATURE_ID} column contains {.strong unique} identifiers!"
            )
          }

        } else {
          private$add_FEATURE_ID()
        }
        cli::cli_alert_success("{.field featureData} is loaded.")
      }

      #-------------------#
      ###   countData   ###
      #-------------------#
      if (!is.null(countData)) {
        private$.countData <- private$check_matrix(countData)
        if (!inherits(private$.countData, "sparseMatrix")) {
          cli::cli_abort(c(
            "Error in {.field countData}:",
            "x" = cli::format_inline("{private$.countData}")
          ))
        }
        cli::cli_alert_success("{.field countData} is loaded.")

        if (is.null(private$.featureData)) {
          private$add_featureData()
          cli::cli_alert_info("Created placeholder {.field featureData}.")
        } else {
          rownames(private$.countData) <- private$.featureData[[ private$.feature_id ]]
        }
      }

      #-------------------#
      ###     sync      ###
      #-------------------#
      private$sync()

      # saves data for reset function
      private$original_data = list(
        counts = private$.countData,
        features = private$.featureData,
        metadata = private$.metaData,
        tree = private$.treeData
      )
    },
    #' @description
    #' This method is very similar to the existing [`clone()`](#method-clone) function, 
    #' except it also creates a new save of the back-up of each class that is created when 
    #' [`new()`](#method-new) is called and can be reverted to via [`reset()`](#method-reset)
    #' 
    #' @param deep A boolean value to create a shallow or deep copy (default: \code{FALSE}).
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #'
    #' obj <- omics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file
    #' )
    #'
    #' # Perform a modification and copy
    #' obj$scale()
    #'
    #' cloned <- obj$copy(deep=TRUE)
    #' cloned$scale(method = "clr")
    #' cloned$reset() # resets to data after clone creation.
    #' 
    #' @return A copy of \code{omics} object.
    copy = function(deep = FALSE) {
      # Base clone
      cloned <- self$clone(deep)

      # Resetting back-up
      cloned$.__enclos_env__$private$original_data <- list(
        counts = private$.countData,
        features = private$.featureData,
        metadata = private$.metaData,
        tree = private$.treeData
      )

      cloned
    },
    #' @description
    #' Displays field parameters of the class via stdout.
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #'
    #' obj <- omics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file
    #' )
    #'
    #' # method 1 to call print function
    #' obj
    #'
    #' # method 2 to call print function
    #' obj$print()
    #'
    #' @return A cli formatted message
    print = function() {
      cli::cli_h3("{.cls {class(self)[1]}} object")
      if (length(private$.metaData) > 0) 
        cli::cli_inform("{.field metaData}: {.val {ncol(private$.metaData)}} variables {cli::symbol$times} {.val {nrow(private$.metaData)}} samples")
      if (length(private$.countData) > 0) 
        cli::cli_inform("{.field countData}: {.val {ncol(private$.countData)}} samples {cli::symbol$times} {.val {nrow(private$.countData)}} features")
      if (length(private$.featureData) > 0)
        cli::cli_inform("{.field featureData}: {.val {ncol(private$.featureData)-1}} attributes {cli::symbol$times} {.val {nrow(private$.featureData)}} features")
      if (length(private$.treeData) > 0)
        cli::cli_inform("{.field treeData}: {.val {length(private$.treeData$tip.label)}} tips {cli::symbol$times} {.val {private$.treeData$Nnode}} nodes")
    },
    #' @description
    #' Upon creation of a new class object (i.e. \code{omics}, \code{metagenomics} or \code{proteomics}) a backup of the original data is created.
    #' Since modification of the object is done by reference and duplicates are not made, it is possible to revert back to the unmodified data via \code{$reset()}.
    #' The methods from the super class \link{omics} also contains a private method to prevent any changes to the original object when using methods,
    #' such as \code{ordination}, \code{alpha_diversity} or \code{foldchange}.
    #' @examples
    #' library(ggplot2)
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' taxa <- omics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file
    #' )
    #'
    #' # Performs modifications
    #' taxa$scale(transform = log2)
    #'
    #' # resets
    #' taxa$reset()
    #'
    #' # An inbuilt reset function prevents unwanted modification to the taxa object.
    #' taxa$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
    #'
    #' @return The object is modified in place and returned (invisibly). If a copy is required, please use [`copy()`](#method-copy) beforehand with \code{obj2 <- obj$copy()}.
    reset = function() {
      if (!is.null(private$original_data)) {
        private$.countData = private$original_data$counts
        private$.featureData = private$original_data$features
        private$.metaData = private$original_data$metadata
        private$.treeData = private$original_data$tree
        invisible(self)
      } else cli::cli_alert_warning("There is no back-up of the data made. This typically happens when the class is not initialized via the {fun. new}.")
    },
    #' @description
    #' Subset `metaData` via the \link[stats]{complete.cases} given a single or vector of column names, then synchronize the object.
    #' @param column A character vector of a single or multiple column names to remove NAs.
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #' 
    #' obj$removeNAs(column = "treatment")
    #' 
    #' @return The object is modified in place and returned (invisibly). If a copy is required, please use [`copy()`](#method-copy) beforehand with \code{obj2 <- obj$copy()}.
    removeNAs = function(column) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (all(is.wholenumber(column)) && length(column) <= length(colnames(private$.metaData)))
        column <- colnames(private$.metaData[column])

      if (!is.character(column))
        cli::cli_abort("{.val column} needs to be a character or an integer.")

      if (!column_exists(column, private$.metaData))
        cli::cli_abort("{.val {column}} does not exist in the {.field metaData} or one of the specified columns is completely empty!")

      ## MAIN
      #--------------------------------------------------------------------#
      private$.metaData <- private$.metaData[stats::complete.cases(private$.metaData[, .SD, .SDcols = column])]
      private$sync()
      invisible(self)
    },
    #' @description
    #' Subset features using a logical expression evaluated in `featureData`, then synchronize the object.
    #' @param expr A logical expression evaluated in the context of `featureData`. Rows for which the expression evaluates to \code{TRUE} are retained.
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #' 
    #' obj$feature_subset(Genus == "Pseudomonas")
    #' 
    #' @return The object is modified in place and returned (invisibly). If a copy is required, please use [`copy()`](#method-copy) beforehand with \code{obj2 <- obj$copy()}.
    feature_subset = function(expr) {
      expr <- base::substitute(expr)

      # Replace all NAs by empty string
      features <- data.table::copy(private$.featureData)
      features[, names(features) := lapply(.SD, function(x) {
        if (is.character(x)) ifelse(is.na(x), "", x) else x
      })]

      ## Evaluate expression in local environment
      rows_to_keep <- tryCatch(
        eval(expr, features, parent.frame()),
        error = function(e) {
          stop(
            sprintf("Failed to evaluate subset expression `%s`: %s",
                    deparse(expr), e$message),
            call. = FALSE
          )
        }
      )
      # Abort if all rows are `FALSE`
      if (base::all(rows_to_keep == FALSE)) {
        cli::cli_abort("The expression resulted in all rows in {.field featureData} to be `FALSE`.")
      }
      ## subset by `rows_to_keep`
      private$.featureData <- private$.featureData[rows_to_keep, , drop = FALSE]
      private$.countData <- private$.countData[rows_to_keep, , drop = FALSE]

      private$sync()
      self$print()
      invisible(self)
    },
    #' @description
    #' Subset samples using a logical expression evaluated in `metaData`, then synchronize the object.
    #' @param expr A logical expression evaluated in the context of `metaData`. Rows for which the expression evaluates to \code{TRUE} are retained.
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #' 
    #' obj$sample_subset(treatment == "tumor")
    #'
    #' @return The object is modified in place and returned (invisibly). If a copy is required, please use [`copy()`](#method-copy) beforehand with \code{obj2 <- obj$copy()}.
    sample_subset = function(expr) {
      expr <- base::substitute(expr)
      metadata <- data.table::copy(private$.metaData)

      ## Evaluate expression in local environment
      rows_to_keep <- tryCatch(
        eval(expr, metadata, parent.frame()),
        error = function(e) {
          stop(
            sprintf("Failed to evaluate subset expression `%s`: %s",
                    deparse(expr), e$message),
            call. = FALSE
          )
        }
      )
      # Abort if all rows are `FALSE`
      if (base::all(rows_to_keep == FALSE)) {
        cli::cli_abort("The expression resulted in all rows in {.field metaData} to be `FALSE`.")
      }

      ## subset by `rows_to_keep`
      private$.metaData <- private$.metaData[rows_to_keep, , drop = FALSE]
      private$.countData <- private$.countData[, private$.metaData[[ private$.sample_id ]], drop = FALSE]

      private$sync()
      self$print()
      invisible(self)
    },
    #' @description
    #' Subset samplepairs by choosing the number of unique pairs in `metaData`, then synchronize the object.
    #' @param num_unique_pairs An integer value to define the number of pairs to subset. The default is NULL, 
    #' meaning the maximum number of unique pairs will be used to subset the data. 
    #' Let's say you have three samples for each pair, then the `num_unique_pairs` will be set to 3.
    #' 
    #' @return The object is modified in place and returned (invisibly). If a copy is required, please use [`copy()`](#method-copy) beforehand with \code{obj2 <- obj$copy()}.
    samplepair_subset = function(num_unique_pairs = NULL) {

      ## Error handling
      #--------------------------------------------------------------------#
      if (!column_exists(private$.samplepair_id, private$.metaData))
        cli::cli_abort("{.val {private$.samplepair_id}} doesn't exist in the {.field metaData}.")

      if (!is.null(num_unique_pairs) && !is.wholenumber(num_unique_pairs))
        cli::cli_abort("{.val num_unique_pairs} must contain integers!")

      ## MAIN
      #--------------------------------------------------------------------#

      counts <- private$.metaData[, .(unique_count = data.table::uniqueN(SAMPLE_ID)), by = SAMPLEPAIR_ID]

      if (is.null(num_unique_pairs)) {
        num_unique_pairs <- counts[, max(unique_count)]
        cli::cli_alert_info("{.val num_unique_pairs} is `NULL`, therefore {.val num_unique_pairs} will be set to {.val {max(counts$unique_count)}}.")
      }

      private$.metaData <- private$.metaData[SAMPLEPAIR_ID %in% counts[unique_count == num_unique_pairs, SAMPLEPAIR_ID]]
      private$sync()
      self$print()
      invisible(self)
    },
    #' @description
    #' Agglomerates features by column, then synchronize the object.
    #' @param feature_rank A character (vector) of columns to aggregate from the `featureData`.
    #' @param feature_filter A character (vector) of characters to remove features via regex substring pattern matching (default: \code{NULL}).
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #' 
    #' obj$feature_merge(feature_rank = c("Kingdom", "Phylum"))
    #' obj$feature_merge(feature_rank = "Genus", feature_filter = c("uncultured", "metagenome"))
    #'
    #' @return The object is modified in place and returned (invisibly). If a copy is required, please use [`copy()`](#method-copy) beforehand with \code{obj2 <- obj$copy()}.
    feature_merge = function(feature_rank, feature_filter = NULL) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank))
        cli::cli_abort("{.val feature_rank} needs to be a character or vector containing characters")

      if (!column_exists(feature_rank, private$.featureData))
        cli::cli_abort("{.val {feature_rank}} does not exist in {.field featureData}!")

      if (!is.null(feature_filter) && !is.character(feature_filter))
        cli::cli_abort("{.val feature_filter} needs to be a character or vector containing characters")

      ## MAIN
      #--------------------------------------------------------------------#
      # Create temporary copies
      counts <- private$.countData
      features <- data.table::copy(private$.featureData[!is.na(private$.featureData[[ feature_rank[1] ]])])

      # set keys
      data.table::setkey(features, FEATURE_ID)

      # Create groups by ID
      grouped_ids <- features[, .(IDs = list(FEATURE_ID)), by = feature_rank]

      ## Skip over singletons
      groups <- base::which(base::lengths(grouped_ids$IDs) > 1)
      grouped_ids$ID_first <- sapply(grouped_ids$IDs, `[[`, 1)

      ## Flatten sparseMatrix to only specific feature ids
      counts <- counts[grouped_ids$ID_first, ]

      ## Sum over multiples 
      for (i in groups) {
        i_id <- grouped_ids$ID_first[i]
        i_group <- grouped_ids$IDs[[i]]
        counts[i_id, ] <- Matrix::colSums(private$.countData[i_group, ])
      }

      # Prepare final self-components
      private$.featureData <- base::unique(features, by = feature_rank)
      
      # Add new agglomerated files
      private$.featureData <- private$.featureData[ 
        base::order(base::match(
          x = private$.featureData[[ private$.feature_id ]], 
          table = grouped_ids$ID_first
        )) 
      ]
      private$.countData <- counts

      # Replaces strings matching feature_filter with NAs
      if (!is.null(feature_filter)) {
        regex_pattern <- paste(feature_filter, collapse = "|")
        for (col in feature_rank) {
          private$.featureData[
            grepl(regex_pattern, get(col), ignore.case = TRUE),
            (col) := NA_character_
          ]
        }
      }

      # Clean up featureData
      empty_strings <- !is.na(private$.featureData[[ feature_rank[1] ]])
      private$.featureData <- private$.featureData[empty_strings, ]
      private$.countData <- private$.countData[empty_strings, ]
      rownames(private$.countData) <- private$.featureData[[ private$.feature_id ]]

      private$sync()
      self$print()
      invisible(self)
    },
    #' @description
    #' Feature scaling on the `countData`. The \code{$scale} function is able to apply transformations element-wise on the positive values, (optional: add pseudocounts) and perform normalisation or standardisation methods.
    #' @param method A character to choose a standardisation/normalisation method, options: `"tss"`, `"clr"`, `"binary"`, `"hellinger"`, `"none"` (default: \code{"tss"}).
    #' \describe{
    #'  \item{\code{"tss"}}{
    #'    Total-sum scaling. Each column is divided by their column sum, converting
    #'    absolute abundance to relative abundances in range from 0 to 1.
    #'  }
    #'  \item{\code{"clr"}}{
    #'    Centered log-ratio-like transformation. Non-zero entries are
    #'    log-transformed using \code{log} and \code{base}, then each feature row is centered by 
    #'    subtracting by their respective geometric mean. 
    #'    Input counts must be strictly positive for his transformation. 
    #'    If \code{pseudocount} is used, then output is equal to \code{vegan::decostand(method = "clr")},
    #'    otherwise the current \code{"clr"} approximates the \code{vegan::decostand(method = "rclr")}, 
    #'    the main difference is that vegan replaces all zero's by NAs and then excludes them when the
    #'    arithmetic mean is taken, whereas zero's in \code{"clr"} are only not used for log-transformation 
    #'    but preserved when taking the arithmetic mean of each feature row.
    #'  }
    #'  \item{\code{"binary"}}{
    #'    Presence-absence transformation. Every stored non-zero count is
    #'    replaced by one, while implicit zeros in the sparse matrix remain zero.
    #'  }
    #'  \item{\code{"hellinger"}}{
    #'    Hellinger transformation. Counts are first total-sum scaled with \code{"tss"}
    #'    to relative abundances, after which the square root of each non-zero value is taken.
    #'  }
    #'  \item{\code{"none"}}{
    #'    No transformation. Typically applied when \code{transform} is used.
    #'  }
    #' }
    #' @param transform A function to apply on the positive values of `countData`, skip standardisation/normalisation with \code{method = "none"} (default: \code{NULL}).
    #' @param base Input for \link[base]{log} to use natural logarithmic scale, log2, log10 or other (default: \code{exp(1)}) in CLR.
    #' @param pseudocount A numeric value to replace zero's (default: \code{NULL}).
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #' # standard relative abundance computation
    #' obj$scale()
    #' 
    #' # CLR
    #' obj$reset()
    #' obj$scale(method = "clr")
    #' 
    #' # transform
    #' obj$reset()
    #' obj$scale(method = "none", transform = log2)
    #' 
    #' @return The object is modified in place and returned (invisibly). If a copy is required, please use [`copy()`](#method-copy) beforehand with \code{obj2 <- obj$copy()}.
    scale = function(method = "tss", transform = NULL, base = exp(1), pseudocount = NULL) {

      ## Nested Functions
      #--------------------------------------------------------------------#
      tss <- function(x) {
        x@x <- x@x / rep(Matrix::colSums(x), base::diff(x@p))
        x
      }

      ## Error handling
      #--------------------------------------------------------------------#
      OPTIONS <- c("tss", "clr", "binary", "hellinger", "none")
      if (!is.character(method) || length(method) != 1) {
        cli::cli_abort("{.val method} needs to contain characters with length of 1.")
      } else if (!method %in% OPTIONS) {
        cli::cli_abort("{.val {method}} is not a valid method. Valid options: <{.val {OPTIONS}}>")
      }

      if (!is.null(pseudocount)) {
        if (!is.numeric(pseudocount) || length(pseudocount) != 1) {
          cli::cli_abort("{.val pseudocount} needs to be a {.cls numeric} type with length of 1.")
        }
      }
        
      if (!is.numeric(base) || length(base) != 1)
        cli::cli_abort("{.val base} needs to be a {.cls numeric} type with length of 1.")

      if (!is.null(transform) && !is.function(transform))
        cli::cli_abort("{.val transform} must be a function!")

      ## MAIN
      #--------------------------------------------------------------------#
      if (!is.null(pseudocount)) 
        private$.countData <- methods::as(private$.countData + pseudocount, "CsparseMatrix")
      
      if (is.function(transform))
        private$.countData@x <- transform(private$.countData@x)

      private$.countData <- switch(
        method,
        "tss" = tss(private$.countData),
        "clr" = {
          ref <- private$.countData
          ref@x <- log(ref@x, base=base)
          ref - Matrix::rowMeans(ref)
        },
        "binary" = {
          ref <- private$.countData
          ref@x[] <- 1
          ref
          },
        "hellinger" = {
          ref <- tss(private$.countData)
          ref@x <- sqrt(ref@x)
          ref
        },
        "none" = private$.countData
      )

      invisible(self)
    },
    #' @description
    #' Rank statistics based on `featureData`. 
    #' Counts the number of (unique) features identified for each column of interest from the `featureData` (excluding the `FEATURE_ID`).
    #' @param feature_ranks A vector of characters that match the `featureData`.
    #' @param unique A boolean value to display only unique entries in `feature_ranks` (default: \code{FALSE}).
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #' 
    #' plt <- obj$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
    #' plt
    #' @return A \link[ggplot2]{ggplot} object.
    rankstat = function(feature_ranks, unique = FALSE) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_ranks))
        cli::cli_abort("{.val feature_ranks} needs to contains characters.")

      if (!column_exists(feature_ranks, private$.featureData))
        cli::cli_abort("Specified {.val {feature_ranks}} do not exist in the {.field featureData}.")

      if (length(feature_ranks) > length(colnames(private$.featureData)))
        feature_ranks <- colnames(private$.featureData[feature_ranks])

      if (!is.logical(unique))
        cli::cli_abort("{.val unique} needs to be either `TRUE` or `FALSE`.")

      ## MAIN
      #--------------------------------------------------------------------#

      # Counts number of ASVs without empty values
      if (unique) {
        values <- private$.featureData[, 
          lapply(.SD, data.table::uniqueN)
          ][, .SD, .SDcols = feature_ranks] 
      } else {
        values <- private$.featureData[, 
          lapply(.SD, function(x) sum(!is.na(x) & x != ""))
          ][, .SD, .SDcols = feature_ranks]  
      }
      
      # Pivot into long table
      long_values <- data.table::melt(
        data = values,
        measure.vars = names(values),
        variable.name = "variable",
        value.name = "counts"
      )

      # Sets order level of taxonomic ranks
      long_values[, variable := factor(variable, levels = base::rev(feature_ranks))]

      # Returns rankstat plot
      return(ggplot2::ggplot(
        data = long_values,
        mapping = ggplot2::aes(
          x = variable,
          y = counts
        )
      ) +
      ggplot2::geom_col(
        fill = "grey",
        colour = "grey15",
        linewidth = 0.25
      ) +
      ggplot2::coord_flip() +
      ggplot2::geom_text(
        mapping = ggplot2::aes(label = counts),
        hjust = -0.1,
        fontface = "bold"
      ) +
      ggplot2::ylim(0, max(long_values$counts)*1.10) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Rank",
        y = "Number of features classified")
      )
    },
    #' @description
    #' Alpha diversity computation, statistical testing and visualisation.
    #' 
    #' It computes the alpha diversity via \link{diversity}, performs a pairwise wilcoxon statistical test given
    #' `groups` via the \link{boxjitter_test}.
    #' 
    #' @param col_name `r lifecycle::badge("deprecated")` This argument has been renamed to `groups` for more clarity.
    #' @param groups A column name in `metaData` to define the contrast of interest for statistical testing and visualisation.
    #' @param metric An alpha diversity metric as input to \link{diversity}, options: `"shannon"`, `"simpson"`, `"invsimpson"` (default: \code{"shannon"}).
    #' @param group_by `r lifecycle::badge("deprecated")` This argument has been renamed to `split_by` for more clarity
    #' @param split_by A column name to split the groups into chunks for grouped statistical test in \link{boxjitter_test} (default: \code{NULL}).
    #' @param Brewer.palID A character name for the palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link{colormap} (default: \code{"Set2"}).
    #' @param evenness A boolean wether to divide diversity by number of species, see \link[vegan]{specnumber} (default: \code{FALSE}).
    #' @param paired A boolean value to perform paired analysis in \link[stats]{wilcox.test} and samplepair subsetting via [`samplepair_subset()`](#method-samplepair_subset)
    #' @param p.adjust.method A character variable to specify the p.adjust.method to be used (default: \code{"fdr"}).
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #' 
    #' plt <- obj$alpha_diversity(groups = "treatment",
    #'                            metric = "shannon")
    #'
    #' @return A list of components: \describe{
    #'  \item{data}{A \link[base]{data.frame} from \link{diversity}.}
    #'  \item{stats}{A pairwise statistics from \link[rstatix]{pairwise_wilcox_test}.}
    #'  \item{plot}{A \link[ggplot2]{ggplot} object.}
    #' }
    alpha_diversity = function(
      col_name = lifecycle::deprecated(),
      groups,
      metric = "shannon",
      Brewer.palID = "Set2",
      group_by = lifecycle::deprecated(),
      split_by = NULL,
      evenness = FALSE,
      paired = FALSE,
      p.adjust.method = "fdr"
    ) {

      ## Lifecycle warning
      if (lifecycle::is_present(col_name)) {
        lifecycle::deprecate_warn(
          when = "1.7.0",
          what = "alpha_diversity(col_name)",
          with = "alpha_diversity(groups)",
          details = "Use `groups` to specify the grouping variable."
        )
        if (missing(groups)) {
          groups <- col_name
        }
      }

      if (lifecycle::is_present(group_by)) {
        lifecycle::deprecate_warn(
          when = "1.7.0",
          what = "alpha_diversity(group_by)",
          with = "alpha_diversity(split_by)",
          details = "Use `split_by` to specify the stratification variable."
        )
        if (missing(split_by)) {
          split_by <- group_by
        }
      }

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(groups) || length(groups) != 1) {
        cli::cli_abort("{.val groups} must be a character and of length 1")
      } else if (!column_exists(groups, private$.metaData)) {
        cli::cli_abort("The specified {.val {groups}} does not exist in the {.field metaData}.")
      }

      if (!is.null(split_by)) {
        if (!is.character(split_by) || length(split_by) != 1) {
          cli::cli_abort("{.val split_by} must be a character and of length 1")
        } else if (!column_exists(split_by, private$.metaData)) {
          cli::cli_abort("The specified {.val {split_by}} does not exist in the {.field metaData}.")
        }
        combined_cols <- c(groups, split_by)
      } else combined_cols <- groups

      if (!is.logical(evenness))
        cli::cli_abort("{.val evenness} can only be a `TRUE` or `FALSE`.")
      
      if (!is.logical(paired))
        cli::cli_abort("{.val paired} can only be a `TRUE` or `FALSE`.")

      if (!c(p.adjust.method %in% stats::p.adjust.methods))
        cli::cli_abort("Specified {.val {p.adjust.method}} is not valid. \nValid options: {.val {p.adjust.methods}}")

      ## MAIN
      #--------------------------------------------------------------------#

      # OUTPUT: Plot list
      plot_list <- list()

      # Save omics class components
      .countData <- private$.countData
      .featureData <- private$.featureData
      .metaData <- private$.metaData
      .treeData <- private$.treeData

      # restore on error
      on.exit({
        private$.countData <- .countData
        private$.featureData <- .featureData
        private$.metaData <- .metaData
        private$.treeData <- .treeData
      }, add = TRUE)

      # Remove NAs from `groups`
      self$removeNAs(groups)

      # Subset by samplepair completion
      if ( paired && column_exists(private$.samplepair_id, private$.metaData) )
        self$samplepair_subset()

      # Alpha diversity based on 'metric'
      div <- data.table::data.table(diversity(x = private$.countData, metric=metric))
      data.table::setnames(x = div, old = "V1", new = metric)
      div[, (combined_cols) := private$.metaData[, .SD, .SDcols = c(combined_cols)]]
      # Adjusts for evenness
      if (evenness) div[[ metric ]] <- div[[ metric ]] / log(vegan::specnumber(div[[ metric ]]))

      # get colors
      colors <- colormap(data = private$.metaData, groups = groups, Brewer.palID = Brewer.palID)

      # Create and saves plots
      plot_list$data <- div
      diversity_plt <- boxjitter_test(
        data = stats::na.omit(div),
        values = metric,
        groups = groups,
        split_by = split_by,
        palette = colors,
        test = "wilcox",
        paired = paired,
        p.adjust.method = p.adjust.method
        )

      plot_list$stats <- as.data.frame(diversity_plt$stats)
      plot_list$plot <- diversity_plt$plot

      return(plot_list)
    },
    #' @description
    #' Creates a long `data.table` of `SAMPLE_ID`s, the top 10-15 most abundant bacteria based on the \code{feature_rank} with
    #' any additional columns specified via \code{add_cols}. The \code{$composition()} method works with aggregated 
    #' values by running first [`feature_merge()`](#method-feature_merge), then sorting for the most abundant bacteria based on
    #' the row sums. The output of \code{$composition()} can be used as input to \link{composition_plot}.
    #' @param feature_rank A column name in `featureData` to aggregate via [`feature_merge()`](#method-feature_merge).
    #' @param feature_filter A character or vector of characters to removes features by regex substring pattern matching via [`feature_merge()`](#method-feature_merge) (default: \code{NULL}).
    #' @param col_name `r lifecycle::badge("deprecated")` This argument has been renamed to `add_cols` for more clarity.
    #' @param add_cols Optional, a column name or vector of multiple column names from `metaData` to add to the `data` list name (default": \code{NULL}).
    #' @param feature_top A wholenumber of the top features to visualize, the max is 15, due to a limit of palettes (default: \code{10}).
    #' @param Brewer.palID A character name for the palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link{colormap} (default: \code{"RdYlBu"}).
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #'
    #' result <- obj$composition(feature_rank = "Genus",
    #'                           feature_filter = c("uncultured"),
    #'                           feature_top = 10)
    #'
    #' plt <- composition_plot(data = result$data,
    #'                         palette = result$palette,
    #'                         feature_rank = "Genus")
    #'
    #' @return A list of components: \describe{
    #'  \item{data}{A \link[data.table]{data.table} of feature compositions.}
    #'  \item{palette}{A \link[stats]{setNames} palette from \link{colormap} matching the top features.}
    #' }
    composition = function(
      feature_rank,
      feature_filter = NULL,
      col_name = lifecycle::deprecated(),
      add_cols = NULL,
      feature_top = 10,
      Brewer.palID = "RdYlBu"
    ) {

      # Handle group_by → split_by
      if (lifecycle::is_present(col_name)) {
        lifecycle::deprecate_warn(
          when = "1.7.0",
          what = "composition(col_name)",
          with = "composition(add_cols)",
          details = "Use `add_cols` to add column names from {.field metaData}."
        )
        if (missing(add_cols)) {
          add_cols <- col_name
        }
      }

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.null(add_cols)) {
        if (!is.character(add_cols)) {
          cli::cli_abort("{.val add_cols} must be a character.")
        } else if (!column_exists(add_cols, private$.metaData)) {
          cli::cli_abort("The specified {.val {add_cols}} does not exist in the {.field metaData}.")
        }
      }

      if (length(feature_top) != 1) {
        cli::cli_abort("{.val feature_top} must be a single element.}")
      } else if (!is.wholenumber(feature_top)) {
        cli::cli_abort("{.val feature_top} must be a whole number.")
      }
      
      if (feature_top > 15) {
        cli::cli_abort("The {.val feature_top} cannot be higher than 15.\n This may lead that colors are difficult to be distinguished for color-blind people, therefore the limit is set to 15.")
      }

      ## MAIN
      #--------------------------------------------------------------------#
      # Copies object to prevent modification of omics class components
      .countData <- private$.countData
      .featureData <- private$.featureData
      .metaData <- private$.metaData
      .treeData <- private$.treeData

      # restore on error
      on.exit({
        private$.countData <- .countData
        private$.featureData <- .featureData
        private$.metaData <- .metaData
        private$.treeData <- .treeData
      }, add = TRUE)

      # Agglomerate by feature_rank
      self$feature_merge(feature_rank = feature_rank, feature_filter = feature_filter)

      # Remove NAs when add_cols is specified
      if (!is.null(add_cols) && length(add_cols) == 1)
        self$removeNAs(add_cols)

      # Converts matrix to data.table
      counts <- matrix_to_dtable(private$.countData)

      # Fetch unfiltered and filtered features
      dt <- counts[, (feature_rank) := private$.featureData[[feature_rank]]]

      # Create row_sums
      dt[, row_sum := rowSums(.SD), .SDcols = !c(feature_rank)]

      # Orders by row_sum in descending order
      data.table::setorder(dt, -row_sum)

      # Subset taxa for visualization
      final_dt <- rbind(dt[1:feature_top][, .SD, .SDcols = !c("row_sum")],
                        dt[(feature_top+1):nrow(dt)][, lapply(.SD, function(x) sum(x)),
                                                                 .SDcols = !c(feature_rank, "row_sum")],
                        fill = TRUE)
      
      # Creates palette
      df_taxa_len <- length(final_dt[[feature_rank]])
      if (df_taxa_len-1 <= 15 && df_taxa_len-1 > 10) {
        chosen_palette <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                            "#920000","#924900","#db6d00","#24ff24","#ffff6d")[1:df_taxa_len-1]
      } else {
        chosen_palette <- RColorBrewer::brewer.pal(df_taxa_len-1, Brewer.palID)
      }
      

      # Add 'Others'
      if (df_taxa_len == feature_top+1) {
        final_dt[nrow(final_dt), (feature_rank)] <- "Other"
        taxa_colors_ordered <- stats::setNames(c(chosen_palette, "lightgrey"), final_dt[[feature_rank]])
      } else {
        taxa_colors_ordered <- stats::setNames(chosen_palette, final_dt[[feature_rank]])
      }

      # Pivoting in long table and factoring feature ranke
      final_long <- data.table::melt(final_dt,
                                     id.vars = c(feature_rank),
                                     variable.factor = FALSE,
                                     value.factor = TRUE)
      # Rename colnames for merge step
      colnames(final_long) <- c(feature_rank, private$.sample_id, "value")

      # Adds metadata columns by user input
      if (!is.null(add_cols)) {
        composition_merged <- base::merge(
          x = final_long,
          y = private$.metaData[, .SD, .SDcols = c(private$.sample_id, add_cols)],
          by = private$.sample_id,
          all = TRUE,
          allow.cartesian = TRUE)
        composition_final <- base::unique(composition_merged)
      } else {
        composition_final <- final_long
      }

      # Factors the melted data.table by the original order of Taxa
      # Important for scale_fill_manual taxa order
      composition_final[[feature_rank]] <- factor(composition_final[[feature_rank]], levels = final_dt[[feature_rank]])

      # returns results as list
      return(
        list(
          data = composition_final,
          palette = taxa_colors_ordered
        )
      )
    },
    #' @description
    #' Compute a distance metric from `countData`
    #' @param metric A dissimilarity metric to be applied on the `countData`, 
    #' options: \describe{
    #'  \item{\code{"bray"}}{\link{bray}}
    #'  \item{\code{"jaccard"}}{\link{jaccard}}
    #'  \item{\code{"cosine"}}{\link{cosine}}
    #'  \item{\code{"manhattan"}}{\link{manhattan}}
    #'  \item{\code{"euclidean"}}{\link{euclidean}}
    #'  \item{\code{"jsd"}}{\link{jsd}}
    #'  \item{\code{"canberra"}}{\link{canberra}}
    #'  \item{\code{"unifrac"}}{\link{unifrac}}
    #'  \item{\code{"aitchison"}}{
    #'    The `"aitchison"` metric uses the \code{$scale(method = "clr")} without any pseudocounts, see [`scale()`](#method-scale),
    #'    followed by euclidean dissimilarity computation.
    #'    
    #'    Ref: Aitchison, J. (1986) The Statistical Analysis of Compositional Data. Chapman and Hall, London, 416 p.
    #'  }
    #' }
    #' @param weighted A boolean value, to use abundances (\code{weighted = TRUE}) or absence/presence (\code{weighted=FALSE}) (default: \code{TRUE}).
    #' @param normalize A boolean value, whether to normalize weighted UniFrac distances to be between 0 and 1. Unweighted UniFrac is always normalized (default: \code{TRUE}).
    #' @param base Input for \link[base]{log} to use natural logarithmic scale, log2, log10 or other (default: \code{exp(1)}).
    #' @param threads A wholenumber, indicating the number of threads to use (default: \code{1}).
    #' @return A column x column \link[stats]{dist} object.
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'     metaData = metadata_file,
    #'     countData = counts_file,
    #'     featureData = features_file
    #' )
    #'
    #' obj$feature_subset(Kingdom == "Bacteria")
    #' dist <- obj$distance(metric = "bray")
    distance = function(metric, weighted = TRUE, threads = 1, normalize = TRUE, base = exp(1)) {

      ## Error handling
      #--------------------------------------------------------------------#
      OPTIONS <- c(
        "bray", "jaccard", "cosine", "manhattan",
        "jsd", "canberra", "unifrac", "euclidean", 
        "aitchison"
        )

      if (missing(metric))
        cli::cli_abort("{.val metric} must be specified!")

      if (!is.character(metric) || length(metric) != 1) {
        cli::cli_abort("{.val metric} needs to be a character with a length of 1")
      } else if (!metric %in% OPTIONS) {
        cli::cli_abort("{.val {metric}} is not a valid {.val metric}. \nValid options: {.val {OPTIONS}}")
      } else if (is.null(private$.treeData) && metric == "unifrac") {
        cli::cli_abort("The specified {.val metric} is invalid since no {.field treeData} is supplied.")
      }

      if (!is.wholenumber(threads))
        cli::cli_abort("{.val threads} need to be a whole number!")        

      ## MAIN
      #--------------------------------------------------------------------#

      # Copies object to prevent modification of omics class components
      .countData <- private$.countData
      .featureData <- private$.featureData
      .metaData <- private$.metaData
      .treeData <- private$.treeData

      # restore on error
      on.exit({
        private$.countData <- .countData
        private$.featureData <- .featureData
        private$.metaData <- .metaData
        private$.treeData <- .treeData
      }, add = TRUE)

      distmat <- switch(
        metric,
        "unifrac" = OmicFlow::unifrac(x = private$.countData, tree = private$.treeData, weighted=weighted, normalize=normalize, threads=threads),
        "manhattan" = OmicFlow::manhattan(x = private$.countData, weighted=weighted, threads=threads),
        "canberra" = OmicFlow::canberra(x = private$.countData, weighted=weighted, threads=threads),
        "jaccard" = OmicFlow::jaccard(x = private$.countData, weighted=weighted, threads=threads),
        "bray" = OmicFlow::bray(x = private$.countData, weighted=weighted, threads=threads),
        "jsd" = OmicFlow::jsd(x = private$.countData, weighted=weighted, threads=threads),
        "cosine" = OmicFlow::cosine(x = private$.countData, weighted=weighted, threads=threads),
        "euclidean" = OmicFlow::euclidean(x = private$.countData, weighted=weighted, threads=threads),
        "aitchison" = {
            self$scale(method = "clr", base=base, pseudocount=NULL)
            OmicFlow::euclidean(x = private$.countData, weighted=weighted, threads=threads)
        }
      )
      return(distmat)
    },
    #' @description
    #' The \code{$ordination()} method either computes a distance matrix by specification of the \code{metric} or 
    #' takes a custom distance matrix via \code{distmat}. It then computes either the pairwise PERMANOVA test on
    #' the defined contrast via \code{group_by} via \link{pairwise_adonis} or the pairwise ANOSIM test via \link{pairwise_anosim}.
    #' 
    #' Which pairwise test to compute is decided via the \code{method} argument, where `"pcoa"` results in
    #' the \link{pairwise_adonis} and `"nmds"` to \link{pairwise_anosim}.
    #' 
    #' It also visualises these results as a scree-plot (only when `"pcoa"` is used), scores-plot (see \link{ordination_plot}) 
    #' and an anova-plot (see \link{plot_pairwise_stats}). These plots, alongisde the distance matrix, principal components
    #' are returned as a list.
    #' 
    #' @param metric A dissimilarity or similarity metric to be applied on the `countData`, see [`distance()`](#method-distance) (default: \code{"bray"}).
    #' @param method Ordination method, options: `"pcoa"` or `"nmds"`, see \link[vegan]{wcmdscale} (default: \code{"pcoa"}).
    #' \describe{
    #'  \item{\code{"pcoa"}}{
    #'    When `"pcoa"` (Principal Coordinate Analysis) is chosen, a pairwise PERMANOVA test via \link{pairwise_adonis} is conducted on 
    #'    the contrast specified by \code{group_by} with n number of permutations, \code{perm}.
    #'  }
    #'  \item{\code{"nmds"}}{
    #'    When `"nmds"` (Non-Metric Multidimensional Scaling) is chosen, a pairwise ANOSIM test via \link{pairwise_anosim} is conducted on 
    #'    the contrast specified by \code{group_by} with n number of permutations, \code{perm}. This method returns 2 dimensions by default,
    #'    and a stress score in the `pcs` list name, that acts as the confidence, see \link[vegan]{metaMDS}. Since only two dimensions
    #'    are returned, there is no `scree_plot` output.
    #'  }
    #' }
    #' @param distmat A custom distance matrix as either \link[stats]{dist} or \link[Matrix]{Matrix} (default: \code{NULL}).
    #' @param group_by A character variable in `metaData` to be used for the \link{pairwise_adonis} or \link{pairwise_anosim} statistical test.
    #' @param weighted A boolean value, whether to compute weighted or unweighted dissimilarities, 
    #' see [`distance()`](#method-distance) (default: \code{TRUE}).
    #' @param normalize A boolean value, wether to normalize weighted UniFrac distances to be between 0 and 1, 
    #' see \link{unifrac} (default: \code{TRUE}).
    #' @param threads A wholenumber, indicating the number of threads to use (default: \code{1}).
    #' @param perm_design A function that takes `metaData` and constructs a permutation design with \link[permute]{how}, 
    #' see [example](https://agusinac.github.io/OmicFlow/articles/visualisation.html#beta-diversity) (default: \code{NULL}).
    #' @param perm A whole number to define the number of permutations in \link[vegan]{adonis2} and \link[vegan]{anosim} (default: \code{999}).
    #' @examples
    #' library("ggplot2")
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #' )
    #'
    #' pcoa_plots <- obj$ordination(metric = "bray",
    #'                              method = "pcoa",
    #'                              group_by = "treatment",
    #'                              weighted = TRUE)
    #' pcoa_plots
    #'
    #' @return A list of components: \describe{
    #'  \item{distmat}{A distance dissimilarity in \link[base]{matrix} format.}
    #'  \item{stats}{A statistical test as a \link[base]{data.frame}.}
    #'  \item{pcs}{principal components as a \link[base]{data.frame}.}
    #'  \item{scree_plot}{A \link[ggplot2]{ggplot} object.}
    #'  \item{anova_plot}{A \link[ggplot2]{ggplot} object.}
    #'  \item{scores_plot}{A \link[ggplot2]{ggplot} object.}
    #' } 
    ordination = function(metric = "bray",
                          method = "pcoa",
                          group_by,
                          distmat = NULL,
                          weighted = TRUE,
                          normalize = TRUE,
                          threads = 1,
                          perm_design = NULL,
                          perm = 999) {

      ## Error handling
      #--------------------------------------------------------------------#
      OPTIONS <- c("pcoa", "nmds")

      if (!is.character(method) || length(method) != 1) {
        cli::cli_abort("{.val method} needs to be a character with a length of 1")
      } else if (!method %in% OPTIONS) {
        cli::cli_abort("{.val {method}} is not a valid {.val method}. \nValid options: {.val {OPTIONS}}")
      }
            
      if (missing(group_by))
        cli::cli_abort("{.val group_by} must be specified!")

      if (!is.character(group_by) || length(group_by) != 1) {
        cli::cli_abort("{.val group_by} needs to be a character with a length of 1")
      } else if (!column_exists(group_by, private$.metaData)) {
        cli::cli_abort("{.val group_by} does not exist in the {.field metaData} or is empty.")
      }

      if (!is.null(perm_design) && !is.function(perm_design))
        cli::cli_abort("{.val perm_design} must be a function.")

      if (!is.wholenumber(perm))
        cli::cli_abort("{.val perm} need to be a whole number.")

      if (!is.null(distmat)) {
        tmp <- distmat
        sample_cols <- private$.metaData[[ private$.sample_id ]]

        if (inherits(tmp, "Matrix") || inherits(tmp, "dist")) {
          if (inherits(tmp, "dist")) {
            tmp <- as.matrix(tmp)
          }

          if (is.null(colnames(tmp))) {
            cli::cli_abort("{.val distmat} doesn't contain any colnames!")
          } else {
            expr <- any(sample_cols %in% colnames(tmp))
            if (!expr) {
              cli::cli_abort("None {.val SAMPLE_ID} from {.field metaData} match the {.val distmat} colnames!")
            }
          }

          rm(tmp)

        } else {
          cli::cli_abort("{.val distmat} need to be {.cls Matrix} or {.cls dist}")
        }
      }

      ## MAIN
      #--------------------------------------------------------------------#
      # Copies object to prevent modification of omics class components
      .countData <- private$.countData
      .featureData <- private$.featureData
      .metaData <- private$.metaData
      .treeData <- private$.treeData

      # restore on error
      on.exit({
        private$.countData <- .countData
        private$.featureData <- .featureData
        private$.metaData <- .metaData
        private$.treeData <- .treeData
      }, add = TRUE)

      # Subset by missing values
      self$removeNAs(group_by)
      if (inherits(distmat, "Matrix")) {
        distmat <- distmat[private$.metaData[[ private$.sample_id ]], private$.metaData[[ private$.sample_id ]]]
        distmat <- as.dist(distmat)
      }

      # Creates a list of plots
      plot_list <- list()

      if (is.null(distmat)) {
        distmat <- self$distance(
          metric = metric,
          weighted = weighted,
          normalize = normalize,
          threads = threads
          )
      }

      plot_list$dist <- as.matrix(distmat)

      # Switch case to compute loading scores
      pcs <- switch(
        method,
        "pcoa" = vegan::wcmdscale(d = distmat,
                                  k = 15,
                                  eig = TRUE),
        "nmds" = vegan::metaMDS(distmat,
                                trace = FALSE,
                                autotransform = FALSE)
      )
      if (!is.null(perm_design)) metadata <- private$.metaData else metadata <- NULL
      # Switch case to compute relevant statistics
      stats_results <- switch(
        method,
        "pcoa" = pairwise_adonis(distmat, groups = private$.metaData[[ group_by ]], perm = perm, perm_design = perm_design, metadata = metadata),
        "nmds" = pairwise_anosim(distmat, groups = private$.metaData[[ group_by ]], perm = perm, perm_design = perm_design, metadata = metadata)
      )
      plot_list$anova_data <- stats_results

      # Data table of loading scores
      df_pcs_points <- data.table::data.table(pcs$points)
      n_samples <- ncol(df_pcs_points)

      if (method == "pcoa") {
        # Normalisation of eigenvalues
        pcs$eig_norm <- unlist(lapply(pcs$eig, function(x) x / sum(pcs$eig) * 100))
        colnames(df_pcs_points) <- paste0("PC", 1:n_samples)

      } else if (method == "nmds") {
        df_pcs_points[, stress := pcs$stress]
      }

      # Adds relevant data
      df_pcs_points[, paste(group_by) := private$.metaData[[ group_by ]] ]
      df_pcs_points[, samples := row.names(df_pcs_points) ]
      plot_list$pcs <- df_pcs_points

      if (method == "pcoa") {
        # Scree plot of first 10 dimensions
        tmp <- data.table::data.table(
          dims = seq(length(pcs$eig_norm[1:n_samples])),
          dims.explained = pcs$eig_norm[1:n_samples]
        )

        plot_list$scree_plot <- ggplot2::ggplot(
          data = tmp,
          mapping = ggplot2::aes(
            x = dims,
            y = dims.explained
            )
          ) +
          ggplot2::geom_col() +
          ggplot2::theme_bw() +
          ggplot2::scale_x_continuous(breaks=seq(1, n_samples, 1)) +
          ggplot2::scale_y_continuous(breaks=seq(0, 100, n_samples)) +
          ggplot2::labs(
            title = paste0("Screeplot of ", length(pcs$eig_norm)," PCs"),
            x = "Principal Components (PCs)",
            y = "dissimilarity explained [%]"
          )

        # PERMANOVA
        plot_list$anova_plot <- plot_pairwise_stats(
          data = stats_results,
          group_col = "pairs",
          stats_col = "F.Model",
          label_col = "p.adj",
          y_axis_title = "Pseudo F test statistic",
          plot_title = "PERMANOVA"
        )
        # Loading score plot
        plot_list$scores_plot <- ordination_plot(
          data = df_pcs_points,
          groups = group_by,
          pair=c("PC1", "PC2"),
          dist_explained = pcs$eig_norm[1:2],
          dist_metric = metric
        )

      } else if (method == "nmds") {
        plot_list$anova_plot <- plot_pairwise_stats(
          data = stats_results,
          group_col = "pairs",
          stats_col = "anosimR",
          label_col = "p.adj",
          y_axis_title = "ANOSIM R statistic",
          plot_title = "ANOSIM"
        )

        plot_list$scores_plot <- ordination_plot(
          data = df_pcs_points,
          groups = group_by,
          pair=c("MDS1", "MDS2"),
          dist_metric = metric
        )
      }

      return(plot_list)
    },
    #' @description
    #' Differential feature expression for both paired and non-paired data.
    #' 
    #' The function performs feature agglomeration via [`feature_merge()`](#method-feature_merge) (which can be disabled via \code{feature_merge=FALSE}), 
    #' subsetting to remove NAs in `condition.group` via [`removeNAs()`](#method-removeNAs) and finding samplepairs when \code{paired = TRUE} is supplied
    #' via [`subset_samplepairs()`](#method-subset_samplepairs).
    #'  
    #' The fold-changes can be computed differently based on the `method` and `aggregate_method` options. 
    #' Transformations of the data should be done beforehand via [`scale()`](#method-scale).
    #' Finally, homogeneity of variance is computed based on the selected `aggregate_method` option. If \code{aggregate_method = "mean"} then the \link[matrixTests]{row_levene} is applied, 
    #' and if \code{aggregate_method = "median"} is used then the \link[matrixTests]{row_brownforsythe} is applied. Any filtering of the results is left to the end-user, 
    #' volcano plots are automatically added to the list name `volcano_plot` via \link{volcano_plot}.
    #' 
    #' @param condition.group A character of an existing column name in `metaData`, wherein the conditions A and B are located.
    #' @param condition_A A character (vector) in `condition.group`.
    #' @param condition_B A character (vector) in `condition.group`.
    #' @param method A character to choose the method of fold-change computation (default: \code{"identity"}). 
    #' \describe{
    #'  \item{\code{"identity"}}{Computes fold-change via \code{log2(condition_A) - log2(condition_B)}, and will handle zero's to prevent `Inf` values. 
    #'  Proportional data is also supported and will be automatically detected.}
    #'  \item{\code{"log"}}{Computes fold-change via \code{condition_A - condition_B}.}
    #' }
    #' How to aggregate the fold-change values for \code{condition_A} and \code{condition_B} can be controlled with the `aggregate_method` argument. 
    #' @param aggregate_method A function to aggregate the fold-change values from multiple samples (default: \code{"median"}). 
    #' options: \describe{
    #'  \item{\code{"mean"}}{
    #'    It takes the arithmetic mean for each feature, separately for `condition_A` and `condition_B`.
    #'  }
    #'  \item{\code{"median"}}{
    #'    It takes the median for each feature, separately for `condition_A` and `condition_B`.
    #'  }
    #'  \item{\code{"geomean"}}{
    #'    It computes the geometric mean of the positive counts, similar to the \link[DESeq2]{estimateSizeFactors()}, when \code{type = "poscount"}.
    #'    The only difference, is that the code is adapted to handle a \link[Matrix]{sparseMatrix}.
    #'  }
    #'  \item{\code{"none"}}{
    #'    No aggreation is applied, the whole matrix is used as it is.
    #'  }  
    #' }
    #' 
    #' @param group_by `r lifecycle::badge("deprecated")` This argument has been renamed to `split_by` for more clarity.
    #' @param split_by A character variable of an existing column in `metaData` to split the table in chunks prior to fold-change computation (default: \code{NULL}). When disabled then column names will end with `_in_all`.
    #' @param feature_merge A boolean value wether to call [`feature_merge()`](#method-feature_merge) (default: \code{FALSE}).
    #' @param feature_rank A column in the `featureData` to use as the feature scope (default: \code{"FEATURE_ID"}).
    #' @param feature_filter A character or vector of characters to remove features via regex substring pattern matching via [`feature_merge()`](#method-feature_merge) (default: \code{NULL}).
    #' @param paired A boolean value, the paired is only applicable when a `SAMPLEPAIR_ID` column exists within the `metaData` (default: \code{FALSE}), see [`samplepair_subset()`](#method-samplepair_subset).
    #' @param pvalue.threshold A numeric value used as a p-value threshold to label and color significant features (default: \code{0.05}), see \link{volcano_plot}.
    #' @param logfold.threshold A numeric value used as a fold-change threshold to label and color significantly expressed features (default: \code{0.06}), see \link{volcano_plot}.
    #' @param abundance.threshold A numeric value used as an abundance threshold to size the scatter dots based on their mean abundance (default: \code{0}), see \link{volcano_plot}.
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #'
    #' obj <- omics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file
    #' )
    #' ## Method 1: using non-log transformation
    #' obj$scale(method = "tss")
    #' 
    #' dfe <- obj$foldchange(
    #'  feature_rank = "Genus",
    #'  feature_merge = TRUE,
    #'  method = "identity",
    #'  aggregate_method = "median",
    #'  condition.group = "treatment",
    #'  condition_A = "tumor",
    #'  condition_B = "healthy"
    #' )
    #' 
    #' ## Method 2: using log transformation
    #' obj$reset()
    #' obj$scale(method = "clr")
    #' 
    #' dfe <- obj$foldchange(
    #'  feature_rank = "Genus",
    #'  feature_merge = TRUE,
    #'  method = "log",
    #'  aggregate_method = "mean",
    #'  condition.group = "treatment",
    #'  condition_A = "tumor",
    #'  condition_B = "healthy"
    #' )
    #' 
    #' @return A list of components: \describe{
    #'  \item{`split_by` subsets}{A \link[base]{matrix} subset from each \code{split_by} separated by \code{condition_A} and \code{condition_B}.}
    #'  \item{data}{A \link[data.table]{data.table} main output of fold-changes between conditions, contains abundance, fold-change, and homogeneity tests.}
    #'  \item{volcano_plot}{A list of \link[ggplot2]{ggplot} plots for each contrast of \code{A vs B}, number of plots depend on \code{split_by} and \code{condition_A} input.}
    #' }
    foldchange = function(
      condition.group,
      condition_A,
      condition_B,
      method = "identity",
      aggregate_method = "median",
      group_by = lifecycle::deprecated(),
      split_by = NULL,
      feature_merge = FALSE,
      feature_rank = "FEATURE_ID",
      feature_filter = NULL,
      paired = FALSE,
      pvalue.threshold = 0.05,
      logfold.threshold = 0.06,
      abundance.threshold = 0
      ) {

      ## Lifecycle warning
      if (lifecycle::is_present(group_by)) {
        lifecycle::deprecate_warn(
          when = "1.7.0",
          what = "foldchange(group_by)",
          with = "foldchange(split_by)"
        )
        if (missing(split_by)) {
          split_by <- group_by
        }
      }

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(condition.group) || length(condition.group) != 1) {
        cli::cli_abort("{.val condition.group} needs to be a character with a length of 1.")
      } else if (!column_exists(condition.group, private$.metaData)) {
        cli::cli_abort("{.val condition.group} does not exist in the {.field metaData} or is empty.")
      }

      if (!is.character(condition_A)) {
        cli::cli_abort("{.val condition_A} needs to be a character.")
      } else if (!any(condition_A %in% private$.metaData[[ condition.group ]])) {
        cli::cli_abort("{.val condition_A} does not exist in {.val condition.group}.")
      }

      if (!is.character(condition_B)) {
        cli::cli_abort("{.val condition_B} needs to be a character.")
      } else if (!any(condition_B %in% private$.metaData[[ condition.group ]])) {
        cli::cli_abort("{.val condition_B} does not exist in {.val condition.group}.")
      }

      OPTIONS <- c("identity", "log")
      if (!is.character(method) || length(method) != 1) {
        cli::cli_abort("{.val method} needs to be a character with a length of 1.")
      } else if (!c(method %in% OPTIONS)) {
        cli::cli_abort("{.val {method}} is not a valid option. \nValid options: {.val {OPTIONS}}")
      }

      OPTIONS <- c("median", "mean", "geomean", "none")
      if (!is.character(aggregate_method) || length(aggregate_method) != 1) {
        cli::cli_abort("{.val aggregate_method} needs to be a character with a length of 1.")
      } else if (!c(aggregate_method %in% OPTIONS)) {
        cli::cli_abort("{.val {aggregate_method}} is not a valid option. \nValid options: {.val {OPTIONS}}")
      }
      
      if (!is.null(split_by)) {
        if (!is.character(split_by) || length(split_by) != 1) {
          cli::cli_abort("{.val split_by} needs to be a character with a length of 1.")
        } else if (!column_exists(split_by, private$.metaData)) {
          cli::cli_abort("The {.val {split_by}} column does not exist in the {.field metaData}.")
        }
      }

      if (!is.logical(feature_merge))
        cli::cli_abort("{.val feature_merge} needs to be either `TRUE` or `FALSE`.")

      if (paired && !column_exists(private$.samplepair_id, private$.metaData)) {
        cli::cli_alert_warning("Paired is set to {.val {paired}} but {.arg SAMPLEPAIR_ID} does not exist in the {.field metaData}.\n Differential feature analysis will continue now with paired set to {.val FALSE}!")
        paired <- FALSE
      }

      ## MAIN
      #--------------------------------------------------------------------#
      ## Adapted from DESeq2 `geoMeanNZ` but re-implemented in sparse format
      geomean <- function(m) {
        ref <- methods::as(m, "CsparseMatrix")
        ref@x <- suppressWarnings(log(ref@x))
        row_means <- exp( Matrix::rowSums(ref, na.rm = TRUE) / ncol(ref) )
        return(row_means)
      }

      # Final output
      output <- list()

      # Copies object to prevent modification of omics class components
      .countData <- private$.countData
      .featureData <- private$.featureData
      .metaData <- private$.metaData
      .treeData <- private$.treeData

      # restore on error
      on.exit({
        private$.countData <- .countData
        private$.featureData <- .featureData
        private$.metaData <- .metaData
        private$.treeData <- .treeData
      }, add = TRUE)

      # Agglomerate features by `feature_rank` and filter unwanted features
      if (feature_merge) {
        self$feature_merge(
          feature_rank = feature_rank,
          feature_filter = feature_filter
        )
      }

      # Subset by missing values
      self$removeNAs(condition.group)

      # Subset by samplepair completion
      if (paired && column_exists(private$.samplepair_id, private$.metaData))
        self$samplepair_subset()

      # Apply `split_by`
      if (!is.null(split_by)) {
        chunks <- base::split(private$.metaData, by = split_by)
      } else {
        chunks <- list(all = private$.metaData)
      }
      group_names <- names(chunks)

      # Create data.tables for results
      foldchange_dt <- data.table::data.table(feature_rank = private$.featureData[[ feature_rank ]])
      colnames(foldchange_dt) <- feature_rank

      ## Add mean abundance
      mat <- as.matrix(private$.countData)
      foldchange_dt[, "median_abun" := matrixStats::rowMedians(mat)]

      for (group_name in group_names) {
        chunk <- chunks[[ group_name ]]

        chunk_mat <- mat[, chunk[[ private$.sample_id ]]]
        condition_labels <- chunk[[ condition.group ]]

        for (i in seq_along(condition_A)) {
          # Subset by `condition_A`
          mat_A <- as.matrix(chunk_mat[, colnames(chunk_mat)[condition_labels %in% condition_A[i]]])
          mat_B <- as.matrix(chunk_mat[, colnames(chunk_mat)[condition_labels %in% condition_B[i]]])

          # save intermediate condition tables
          output[[paste0(group_name, "_", condition_A[i])]] <- mat_A
          output[[paste0(group_name, "_", condition_B[i])]] <- mat_B

          ## Get mean/median for each condition
          ### condition A
          if (ncol(mat_A) == 1) {
            row_A <- mat_A
          } else {
            row_A <- switch(
              aggregate_method, 
              "mean" = matrixStats::rowMeans2(mat_A), "median" = matrixStats::rowMedians(mat_A),
              "geomean" = geomean(mat_A), "none" = mat_A
              )
          }

          ### condition B
          if (ncol(mat_B) == 1) {
              row_B <- mat_B
            } else {
              row_B <- switch(
                aggregate_method, 
                "mean" = matrixStats::rowMeans2(mat_B), "median" = matrixStats::rowMedians(mat_B),
                "geomean" = geomean(mat_B), "none" = mat_B
              )
            }
          
          ## computing fold-change based on method
          if (method == "log") {
            fc_res <- row_A - row_B

          ## Computing fold-change by division
          } else if (method == "identity") {
            max_val <- base::max(row_A)

            if (is.null(dim(row_A))) {
              fc_res <- numeric(length(row_A))
            } else {
              fc_res <- matrix(
                data = 0, 
                nrow = nrow(row_A), 
                ncol = ncol(row_A),
                dimnames = list(rownames(row_A), colnames(row_A))
              )
            }
            # Find zero's to prevent Inf
            both_zero <- row_A == 0 & row_B == 0
            row_A_zero <- row_A == 0 & row_B != 0
            row_B_zero <- row_A != 0 & row_B == 0
            both_non_zero <- row_A != 0 & row_B != 0

            # Compute log2 fold change
            ## In this case we use log2(A) - log2(B)
            fc_res[both_zero] <- 0
            fc_res[row_A_zero] <- row_A[row_A_zero] - log2(row_B[row_A_zero])
            fc_res[row_B_zero] <- log2(row_A[row_B_zero]) - row_B[row_B_zero]
            fc_res[both_non_zero] <- log2(row_A[both_non_zero]) - log2(row_B[both_non_zero])

            # Reverse flipped values with zero's based on max_val
            if (max_val < 1.0) {
              fc_res[row_A_zero] <- fc_res[row_A_zero] * -1
              fc_res[row_B_zero] <- fc_res[row_B_zero] * -1
            }
          }

          ## Combine fold-change with main table
          foldchange_dt <- cbind(foldchange_dt, fc_res)
          
          result_col_title <- paste0(condition_A[i], "_vs_", condition_B[i], "_in_", group_name)
          if (aggregate_method != "none") {
            colnames(foldchange_dt)[grepl("fc_res", colnames(foldchange_dt))] <- paste0("fold-change_", result_col_title)
          }

          ## Computing wilcox test
          if (paired) {
            foldchange_dt[, (paste0("pvalue_wilcox-paired_", result_col_title)) := matrixTests::row_wilcoxon_paired(x = mat_A, y = mat_B)$pvalue]
          } else {
            foldchange_dt[, (paste0("pvalue_wilcox_", result_col_title)) := matrixTests::row_wilcoxon_twosample(x = mat_A, y = mat_B)$pvalue]
          }
          
          ## Compute homogeneity of variance test based on selected `aggregate_method`
          if (aggregate_method == "geomean") aggregate_method <- "mean"

          combined_mat <- cbind(mat_A, mat_B)
          combined_labels <- c(rep(paste0(condition_A), ncol(mat_A)), rep(paste0(condition_B), ncol(mat_B)))
          homogeneity_test <- switch(
            aggregate_method,
            "mean" = matrixTests::row_levene(x = combined_mat, g = combined_labels),
            "median" = matrixTests::row_brownforsythe(x = combined_mat, g = combined_labels)
          )
          if (!is.null(homogeneity_test)) {
            foldchange_dt[, (paste0("homogeneity_test_statistic_", result_col_title)) := homogeneity_test$statistic]
            foldchange_dt[, (paste0("homogeneity_test_pvalue_", result_col_title)) := homogeneity_test$pvalue]
          }
        }
      }
      output$data <- foldchange_dt

      #----------------------#
      # Visualization        #
      #----------------------#

      # Create & save volcano plot if method is not `"none"`
      if (aggregate_method != "none") {
        colnames_dfe <- colnames(foldchange_dt)
        diff_columns <- colnames_dfe[grepl("fold-change", colnames_dfe)]
        pvalue_columns <- colnames_dfe[grepl("pvalue_wilcox", colnames_dfe)]
        abun_column <- colnames_dfe[grepl("abun", colnames_dfe)]
        n_diff_columns <- length(diff_columns)

        output$volcano_plot <- lapply(1:n_diff_columns, function(i) {
          volcano_plot(
            data = foldchange_dt,
            logfold_col = diff_columns[i],
            pvalue_col = pvalue_columns[i],
            feature_rank = feature_rank,
            abundance_col = abun_column,
            pvalue.threshold = pvalue.threshold,
            logfold.threshold = logfold.threshold,
            abundance.threshold = abundance.threshold,
            label_A = condition_A,
            label_B = condition_B
          ) + ggplot2::labs(
              subtitle = paste0(
                "Attribute: ", condition.group,
                ", test: ", ifelse(paired, "Wilcox signed rank test", "Mann-Whitney U test")
                )
            )
        })
      }

      return(output)
    },
    #' @description
    #' Automated Omics Analysis only works if there are column names in the `metaData` with the prefix \code{"CONTRAST_"},
    #' which is detected via sub-string pattern matching. The transformation of the data is left to the user via the
    #' [`scale()`](#method-scale) method. The \code{$autoFlow()} incorporates the following functions: \describe{
    #'  \item{\code{$rankstat()}}{
    #'    This function is only applied when \code{feature_ranks} is specified, useful in \link{metagenomics},
    #'    see [`rankstat()`](#method-rankstat) for more information.
    #'  }
    #'  \item{\code{$alpha_diversity()}}{
    #'    This function only works with the `metagenomics` class, uses `"shannon"`
    #'    metric, and works for both paired and non-paired data, applied on each `"CONTRAST_"` column,
    #'    see [`alpha_diversity()`](#method-alpha_diversity) for more information.
    #'  }
    #'  \item{\code{$ordination()}}{
    #'    This function is applied on all classes, multiple \code{distance_metrics} can be supplied,
    #'    and only the Principal Component method (\code{$ordination(method = "pcoa")}) is used.
    #'    Additional arguments can be controlled with \code{distance_metrics}, \code{distmat},
    #'    \code{weighted}, \code{perm}, \code{threads}, see [`ordination()`](#method-ordination)
    #'    for more information. **NOTE**: \code{perm_design} is not supported in [`autoFlow()`](#method-autoFlow).
    #'  }
    #'  \item{\code{$composition()}}{
    #'    This function is applied on the significant contrasts that are obtained from either
    #'    \code{$alpha_diversity()} or \code{$ordination()}. Multiple \code{feature_contrasts} can be
    #'    supplied, by default it only uses the \code{feature_contrasts = "FEATURE_ID"}. The 
    #'    composition plots are created by default via the \code{composition_plot(group_by)} with
    #'    \code{feature_top = 15}. Additional arguments can be controlled with \code{feature_filter},
    #'    see [`composition()`](#method-composition) for more information.
    #'  }
    #'  \item{\code{$foldchange()}}{
    #'    This function is applied on the significant contrasts that are obtained from either 
    #'    \code{$alpha_diversity()} or \code{$ordination()}. Multiple \code{feature_contrasts} can be
    #'    supplied, by default it only uses the \code{feature_contrasts = "FEATURE_ID"}. Additional arguments
    #'    can be controlled with \code{fc_method}, \code{aggregate_method}, \code{pvalue.threshold}, 
    #'    \code{logfold.threshold}, \code{abundance.threshold}, 
    #'    see [`foldchange()`](#method-foldchange) for more information.
    #'  }
    #' }
    #' The above functions are called in order, only significant contrasts are further used to create
    #' composition and volcano plots. The \code{autoFlow} is intended as a first screening of certain
    #' contrasts and extra analysis, such as the use of a permutation design are not supported and should
    #' be applied by the user.
    #' 
    #' The \code{report} argument controls if an `rmarkdown` document is created, and with \code{filename}
    #' the output filename can be specified.
    #' 
    #' @param feature_contrast A character vector of feature columns in the `featureData` to aggregate via [`feature_merge()`](#method-feature_merge) (default: \code{"FEATURE_ID"}).
    #' @param feature_filter A character vector to filter unwanted features, see [`feature_merge()`](#method-feature_merge) (default: \code{NULL}).
    #' @param feature_ranks A character vector as input to [`rankstat()`](#method-rankstat) (default: \code{NULL}).
    #' @param distance_metrics A character vector specifying what (dis)similarity metrics to use in [`ordination()`](#method-ordination) (default: \code{c("bray")}).
    #' When you are working with log-transformed data it is advised to use the `"euclidean"`.
    #' @param distmat A path to an existing file or a dense/sparse \link[Matrix]{Matrix} format to be used in [`ordination()`](#method-ordination) (default: \code{NULL}).
    #' @param weighted A boolean value, whether to compute weighted or unweighted dissimilarities in [`ordination()`](#method-ordination) (default: \code{TRUE}).
    #' @param fc_method A character to choose the method in [`foldchange()`](#method-foldchange) (default: \code{"identity"}). 
    #' \describe{
    #'  \item{\code{"identity"}}{Computes fold-change via \code{log2(condition_A) - log2(condition_B)}, and will handle zero's to prevent `Inf` values. 
    #'  Proportional data is also supported and will be automatically detected.}
    #'  \item{\code{"log"}}{Computes fold-change via \code{condition_A - condition_B}.}
    #' }
    #' @param aggregate_method A function to aggregate the matrix values in [`foldchange()`](#method-foldchange) by taking e.g. 
    #' the \code{median} of `condition_A` and `condition_B` prior to substraction, or in the case of \code{method = "identity"} 
    #' to take the median of the fold-change \code{log2(median(A)) - log2(median(B))} (default: \code{median}).
    #' @param pvalue.threshold A numeric value used as a p-value threshold to label and color significant features in [`foldchange()`](#method-foldchange) (default: \code{0.05}).
    #' @param logfold.threshold A numeric value used as a fold-change threshold to label and color significantly expressed features in [`foldchange()`](#method-foldchange) (default: \code{0.06}).
    #' @param abundance.threshold A numeric value used as an abundance threshold to size the scatter dots based on their mean abundance in [`foldchange()`](#method-foldchange) (default: \code{0}).
    #' @param perm A wholenumber, number of permutations to compare against the null hypothesis of \link[vegan]{adonis2} or \link[vegan]{anosim} in [`ordination()`](#method-ordination) (default: \code{999}).
    #' @param threads Number of threads to use, only used in [`distance()`](#method-distance) when distmat is not supplied (default: \code{1}).
    #' @param report A boolean value to create a HTML markdown report (default: \code{FALSE}). If \code{FALSE} a nested list of the plots and data is returned.
    #' @param filename A character to name the HTML report to be saved in the current working directory (default: \code{paste0(getwd(), "/report.html")}). The \code{getwd()} is required for rmarkdown to save it in the right path.
    #' 
    #' @return if \code{report = TRUE} a HTML report is generated, otherwise a list containing: \describe{
    #'  \item{plots}{
    #'    - `rankstat_plot` as a \link[ggplot2]{ggplot}.
    #'    - `alpha_div_plots` as a list of \link[ggplot2]{ggplot}.
    #'    - `pcoa_plots` as a nested matrix of lists with \link[ggplot2]{ggplot} (`contrasts` x `distance_metrics`).
    #'    - `composition_plots` as a nested matrix of lists with \link[ggplot2]{ggplot} (`contrasts` x `feature_contrasts`).
    #'    - `Log2FC_plots` as a nested matrix of lists with \link[ggplot2]{ggplot} (`contrasts` x `feature_contrasts`).
    #'  }
    #' 
    #'  \item{data}{
    #'    - `alpha_div_data` as a list of lists of `data` and `stats`.
    #'    - `pcoa_data` as a nested matrix of lists of `stats`, `dist_mat` and `pcs` (`contrasts` x `distance_metrics`).
    #'    - `composition_data` as a nested matrix of lists of `data` (`contrasts` x `feature_contrasts`).
    #'    - `Log2FC_data` as a nested matrix of lists of `data` (`contrasts` x `feature_contrasts`).
    #'  }
    #' } 
    #' 
    autoFlow = function(
      feature_contrast = "FEATURE_ID",
      feature_filter = NULL,
      feature_ranks = NULL,
      distance_metrics = c("bray"),
      distmat = NULL,
      weighted = TRUE,
      fc_method = "log",
      aggregate_method = "median",
      pvalue.threshold = 0.05,
      logfold.threshold = 0.06,
      abundance.threshold = 0,
      perm = 999,
      threads = 1,
      report = TRUE,
      filename = paste0(getwd(), "/report.html")
    ) {
    ## Error handling
    #--------------------------------------------------------------------#

    if (!is.character(feature_contrast)) {
      cli::cli_abort("{.val feature_contrast} needs to contain characters.")
    } else if (!column_exists(feature_contrast, private$.featureData)) {
      cli::cli_abort("{.val {feature_contrast}} does not exist in {.field featureData}!")
    }
    
    if (!is.logical(report))
      cli::cli_abort("{.val report} needs to be either `TRUE` or `FALSE`.")

    if (!is.character(filename) || length(filename) != 1)
      cli::cli_abort("{.val filename} needs to be a character with a length of 1")
      
    ## MAIN
    #--------------------------------------------------------------------#
    is_empty = function(obj) {
      if (length(obj) == 0) {
        return(NULL)
      } else if (length(obj) > 0) {
        keep_cells <- sapply(obj, function(x) is.null(x))
        obj <- obj[!keep_cells]
        if (length(obj) == 0) {
          return(NULL)
        } else return(obj)
      }
    }

    # Creates empty plots and data list
    plots <- list()
    data <- list()
    
    # Save omics class components
    .countData <- private$.countData
    .featureData <- private$.featureData
    .metaData <- private$.metaData
    .treeData <- private$.treeData

    # restore on error
    on.exit({
      private$.countData <- .countData
      private$.featureData <- .featureData
      private$.metaData <- .metaData
      private$.treeData <- .treeData
    }, add = TRUE)

    # Collect columns: CONTRAST_ and VARIABLE_
    metacols <- colnames(private$.metaData)

    CONTRAST_data <- private$.metaData[, .SD, .SDcols = grepl("CONTRAST_", metacols)]
    CONTRAST_names <- colnames(CONTRAST_data)

    VARIABLE_data <- private$.metaData[, .SD, .SDcols = grepl("VARIABLE_", metacols)]
    VARIABLE_names <- colnames(VARIABLE_data)

    if (ncol(CONTRAST_data) == 0)
      cli::cli_abort("No columns with prefix {.val CONTRAST} found.. Did you forgot to add a prefix?")

    #---------------------------------------------#
    # Perform standard visualizations             #
    #---------------------------------------------#
    #
    # CONTRAST
    #
    feature_nrow <- length(feature_contrast)
    CONTRAST_ncol <- length(CONTRAST_data)
    VARIABLE_ncol <- length(VARIABLE_data)

    # Standard rank stats
    if (!is.null(feature_ranks)) {
      plots$rankstat_plot <- self$rankstat(feature_ranks)
    }

    # Main loop
    if (CONTRAST_ncol > 0) {

      # Load custom distance matrix if supplied
      if (!is.null(distmat)) {
        distmat <- private$check_matrix(distmat)
        if (!inherits(distmat, "sparseMatrix")) {
          cli::cli_abort(c(
            "Error in {.val distmat}:",
            "x" = cli::format_inline("{distmat}")
          ))
        }

        distmat <- distmat[private$.metaData[[private$.sample_id]], private$.metaData[[private$.sample_id]]]
      }

      # Initialize plot containers
      composition_plots <- matrix(list(), CONTRAST_ncol, feature_nrow)
      Log2FC_plots <- matrix(list(), CONTRAST_ncol, feature_nrow)
      alpha_div_plots <- list()
      metrics_nrow <- length(distance_metrics)
      pcoa_plots <- matrix(list(), CONTRAST_ncol, metrics_nrow)

      # Initialize data containers
      composition_data <- matrix(list(), CONTRAST_ncol, feature_nrow)
      Log2FC_data <- matrix(list(), CONTRAST_ncol, feature_nrow)
      alpha_div_data <- list()
      pcoa_data <- matrix(list(), CONTRAST_ncol, metrics_nrow)

      for (i in 1:CONTRAST_ncol) {
        col_name <- CONTRAST_names[i]
        conditions <- NULL
        cli::cli_alert_info(paste0("Processing ... column: ", col_name, " \n"))

        #--------------------------------------------------------------------#
        ## Alpha diversity
        #--------------------------------------------------------------------#
        if ("metagenomics" %in% class(self)) {
          res <- tryCatch(
            {
              # Default attempt
              self$alpha_diversity(
                groups = col_name,
                metric = "shannon",
                paired = ifelse(column_exists(private$.samplepair_id, private$.metaData), TRUE, FALSE)
              )
            },
            error = function(e) {
              cli::cli_alert_warning("{.arg alpha_diversity} with {.val paired=TRUE} failed. Retrying with {.val paired=FALSE}.")

            # Retry with paired = FALSE
            res2 <- tryCatch(
              self$alpha_diversity(
                groups = col_name,
                metric = "shannon",
                paired = FALSE
              ),
              error = function(e2) {
                cli::cli_alert_info("Skipping {.arg alpha_diversity}, which failed due to an error: {.val {e2}}.")
                NULL
                }
              )
              res2
            }
          )

          if (!is.null(res)) {
            ## Save plots & data
            alpha_div_plots[[i]] <- res$plot
            alpha_div_data[[i]] <- list(data = res$data, stats = res$stats)
            
            ### Identify significant groups for composition plots & volcano plots
            signif_pairs <- res$stats[res$stats$p.adj < pvalue.threshold, ][c("group1", "group2")]
            if (nrow(signif_pairs) > 0)
              conditions <- signif_pairs
          }
        }

        #--------------------------------------------------------------------#
        ## Beta diversity
        #--------------------------------------------------------------------#

        for (j in 1:metrics_nrow) {
          if (inherits(distmat, "Matrix")) {
            res <- self$ordination(
              distmat = distmat,
              method = "pcoa",
              perm = perm,
              group_by = col_name
              )
          } else {
            res <- self$ordination(
              metric = distance_metrics[j],
              method = "pcoa",
              group_by = col_name,
              weighted = weighted,
              perm = perm,
              threads = threads
              )
          }
          
          ## Save plots and identify significant groups for composition plots & volcano plots
          signif_pairs <- res$anova_data[res$anova_data$p.adj < pvalue.threshold, ]
          if (nrow(signif_pairs) > 0) {
            pairs_split <- strsplit(as.character(signif_pairs$pairs), " vs ")
            
            # Create group1 and group2 columns from split
            signif_pairs$group1 <- sapply(pairs_split, `[`, 1)
            signif_pairs$group2 <- sapply(pairs_split, `[`, 2)
            
            signif_pairs <- signif_pairs[c("group1", "group2")]
            
            conditions <- combine_conditions(conditions, signif_pairs)
          }
          
          ### Store plot and data
          pcoa_plots[[i, j]] <- patchwork::wrap_plots(res[c("scree_plot", "anova_plot", "scores_plot")],
                                                      nrow = 1) +
            patchwork::plot_layout(widths = c(rep(5, 3)),
                                   guides = "collect")
          pcoa_data[[i, j]] <- list(
            stats = res$anova_data,
            dist_mat = res$dist,
            pcs = res$pcs
          )

          # Creates temporary plot results for NMDS
          if (inherits(distmat, "Matrix")) {
            res <- self$ordination(
              distmat = distmat,
              method = "nmds",
              group_by = col_name,
              perm = perm
              )
          } else {
            res <- self$ordination(
              metric = distance_metrics[j],
              method = "nmds",
              group_by = col_name,
              weighted = weighted,
              perm = perm,
              threads = threads
              )
          }

          ## Save plots and identify significant groups for composition plots & volcano plots
          signif_pairs <- res$anova_data[res$anova_data$p.adj < pvalue.threshold, ]
          if (nrow(signif_pairs) > 0) {
            pairs_split <- strsplit(as.character(signif_pairs$pairs), " vs ")
            
            # Create group1 and group2 columns from split
            signif_pairs$group1 <- sapply(pairs_split, `[`, 1)
            signif_pairs$group2 <- sapply(pairs_split, `[`, 2)
            
            signif_pairs <- signif_pairs[c("group1", "group2")]
            
            conditions <- combine_conditions(conditions, signif_pairs)
          }      
        }
      
        #--------------------------------------------------------------------#
        ## Feature composition & FOLDCHANGE
        #--------------------------------------------------------------------#

        for (j in 1:feature_nrow) {
          if (!any(private$.countData@x < 0, na.rm = TRUE)) {
            # Creates composition long table
            res <- self$composition(
              feature_rank = feature_contrast[j],
              feature_filter = feature_filter,
              feature_top = 15,
              col_name = col_name
              )
            # Creates composition ggplot and stores plot with data
            composition_plots[[i, j]] <- composition_plot(
              data = res$data,
              palette = res$palette,
              feature_rank = feature_contrast[j],
              group_by = col_name
              )
            composition_data[[i, j]] <- list(data = res$data)
          } else {
            cli::cli_alert_info("Skipping {.function composition} method due to the detection of negative values.")
          }
          
          if (!is.null(conditions) && nrow(conditions) > 0) {

            dfe <- tryCatch(
              {
              # Default attempt
              self$foldchange(
                feature_merge = ifelse(feature_contrast[j] == "FEATURE_ID", FALSE, TRUE),
                feature_rank = feature_contrast[j],
                feature_filter = feature_filter,
                aggregate_method = aggregate_method,
                method = fc_method,
                paired = ifelse(column_exists(private$.samplepair_id, private$.metaData), TRUE, FALSE),
                condition.group = col_name,
                condition_A = c(conditions$group1),
                condition_B = c(conditions$group2),
                pvalue.threshold = pvalue.threshold,
                abundance.threshold = abundance.threshold,
                logfold.threshold = logfold.threshold
                )
              },
              error = function(e) {
                cli::cli_alert_warning("{.function foldchange} with `paired=TRUE` failed. Retrying with `paired=FALSE`.")
                self$foldchange(
                  feature_merge = ifelse(feature_contrast[j] == "FEATURE_ID", FALSE, TRUE),
                  feature_rank = feature_contrast[j],
                  feature_filter = feature_filter,
                  aggregate_method = aggregate_method,
                  method = fc_method,
                  paired = FALSE,
                  condition.group = col_name,
                  condition_A = c(conditions$group1),
                  condition_B = c(conditions$group2),
                  pvalue.threshold = pvalue.threshold,
                  abundance.threshold = abundance.threshold,
                  logfold.threshold = logfold.threshold
                )
              }
            )
            Log2FC_plots[[i, j]] <- patchwork::wrap_plots(dfe$volcano_plot, nrow=1)
            Log2FC_data[[i, j]] <- list(data = dfe$data)
          }
        }
      }
      
      # Checks if plots aren't empty
      plots$alpha_div_plots <- is_empty(alpha_div_plots)
      plots$composition_plots <- is_empty(composition_plots)
      plots$Log2FC_plots <- is_empty(Log2FC_plots)
      plots$pcoa_plots <- is_empty(pcoa_plots)

      # Checks if data aren't empty
      data$composition_data <- is_empty(composition_data)
      data$Log2FC_data <- is_empty(Log2FC_data)
      data$alpha_div_data <- is_empty(alpha_div_data)
      data$pcoa_data <- is_empty(pcoa_data)
    }
    
    #--------------------------------------------------------------------#
    ## CREATING REPORT
    #--------------------------------------------------------------------#
    if (report) {
      # Locate the template Rmd and CSS within the installed package
      rmd_path <- system.file("report.Rmd", package = "OmicFlow")
      css_path <- system.file("styles.css", package = "OmicFlow")

      ## To bypass R CMD error and define for docker
      knit_dir <- dirname(filename)
      
      rmarkdown::render(
        input = rmd_path,
        output_file = filename,
        intermediates_dir = knit_dir,
        knit_root_dir = knit_dir,
        output_options = list(css = css_path)
      )
    } else {
      return(list(
        plots = plots,
        data = data
      ))
    }
  }
  ),
  private = list(

    # Private data fields
    #-------------------------#
    .countData = NULL,
    .featureData = NULL,
    .metaData = NULL,
    .treeData = NULL,
    .valid_schema = NULL,
    .feature_id = "FEATURE_ID",
    .sample_id = "SAMPLE_ID",
    .samplepair_id = "SAMPLEPAIR_ID",
    original_data = list(),

    # Validate metadata schema
    #---------------------------------------------------------#
    validate = function() {
      # Creates temporary json file from `metaData`
      tmp_json <- base::tempfile(fileext = ".json")

      yyjsonr::write_json_file(
        x = private$.metaData,
        filename = tmp_json
      )

      # Check against schema
      private$.valid_schema <- jsonvalidate::json_validate(
        tmp_json,
        system.file("metadata_schema.json", package = "OmicFlow"),
        engine = "ajv",
        verbose = TRUE,
        error = FALSE,
        strict = TRUE
      )

      unlink(tmp_json)
      invisible(self)
    },
    # Function for synchronization of private data fields
    #---------------------------------------------------------#
    sync = function() {
      if (!is.null(private$.metaData)) {
        if (!column_exists(private$.sample_id, private$.metaData))
          return("{private$.sample_id} doesn't exist in {.field metaData}.")

        private$.metaData <- private$.metaData[, lapply(.SD, function(x) ifelse(x == "", NA, x)),
                                        .SDcols = colnames(private$.metaData)]

        colnames(private$.metaData) <- gsub("\\s+", "_", colnames(private$.metaData))  

        # Keep only common samples based on metaData
        if (!is.null(private$.countData)) {
          private$.countData <- private$check_matrix(private$.countData)
          if (!inherits(private$.countData, "sparseMatrix")) {
            cli::cli_abort(c(
              "Error in {.field countData}:",
              "x" = cli::format_inline("{private$.countData}")
            ))
          }
          common_samples <- base::intersect(private$.metaData[[ private$.sample_id ]], colnames(private$.countData))

          if (length(common_samples) == 0)
            cli::cli_abort("None SAMPLE_IDs are matching, check if {.val SAMPLE_ID} are matching the colnames in {.field countData}!")

          private$.countData <- private$.countData[, common_samples, drop = FALSE]
          private$.metaData <- private$.metaData[private$.metaData[[ private$.sample_id ]] %in% common_samples, ]
        }
      }

      if (!is.null(private$.featureData)) {
        if (!column_exists(private$.feature_id, private$.featureData))
          cli::cli_abort("{private$.feature_id} doesn't exist in {.field featureData}.")

        private$.featureData <- private$check_table(private$.featureData)
        if (!data.table::is.data.table(private$.featureData)) {
          cli::cli_abort(c(
            "Error in {.field featureData}:",
            "x" = cli::format_inline("{private$.featureData}")
          ))
        }

        colnames(private$.featureData) <- gsub("\\s+", "_", colnames(private$.featureData))

        # Keep only common tips based on treeData
        if (!is.null(private$.treeData)) {
          common_tips <- base::intersect(private$.treeData$tip.label, private$.featureData[[ private$.feature_id ]])

          if (length(common_tips) == 0)
            cli::cli_abort("None FEATURE_IDs are matching, check if {.val FEATURE_ID} matches the tip labels in {.field treeData}!")

          private$.treeData <- ape::keep.tip(private$.treeData, common_tips)
          private$.featureData <- private$.featureData[private$.featureData[[ private$.feature_id ]] %in% common_tips, ]
        }

        # Keep only common features based on countData
        if (!is.null(private$.countData)) {
          common_features <- base::intersect(private$.featureData[[ private$.feature_id ]], rownames(private$.countData))

          if (length(common_features) == 0)
            cli::cli_abort("None FEATURE_IDs are matching, check if {.val FEATURE_ID} matches the rownames in {.field countData}!")

          private$.featureData <- private$.featureData[private$.featureData[[ private$.feature_id ]] %in% common_features, ]
          private$.countData <- private$.countData[common_features, ]
          private$removeZeros()
        }
      } else if (!is.null(private$.countData)) {
        private$add_featureData()
        cli::cli_alert_info("Placeholder {.field featureData} created.")
      }
    },
    # Helper function to remove empty rows/columns
    #--------------------------------------#
    removeZeros = function() {
      keep_cols <- base::diff(private$.countData@p) > 0
      keep_rows <- base::diff(Matrix::t(private$.countData)@p) > 0

      private$.countData <- private$.countData[keep_rows, keep_cols]
      private$.metaData <- private$.metaData[keep_cols, ]
      private$.featureData <- private$.featureData[keep_rows]

      if (!is.null(private$.treeData))
        private$.treeData <- ape::keep.tip(private$.treeData, private$.featureData[[ private$.feature_id ]])
    },
    # Adds a `featureData`
    #--------------------------------------#
    add_featureData = function() {
      private$.featureData <- data.table::data.table()
      countData_with_rownames <- rownames(private$.countData)

      if (is.null(countData_with_rownames)) {
        private$.featureData[, (private$.feature_id) := paste0("feature_", 1:nrow(private$.countData))]
        rownames(private$.countData) <- private$.featureData[[ private$.feature_id ]]
      } else {
        private$.featureData[, (private$.feature_id) := countData_with_rownames]
      }
      private$add_FEATURE_ID()       
    },
    # Adds or re-aligns `FEATURE_ID` as first column
    #--------------------------------------#
    add_FEATURE_ID = function() {
      if (!column_exists(private$.feature_id, private$.featureData)) {
        FEATURE_ID <- paste0("feature_", 1:nrow(private$.featureData))
        private$.featureData[, (private$.feature_id) := FEATURE_ID]
      }
      data.table::setcolorder(
        x = private$.featureData, 
        neworder = c(private$.feature_id, base::setdiff(private$.feature_id, colnames(private$.featureData)))
      )
    },
    # Extension of `countData`
    # Atm this function is only applied when `countData` is a `data.table`
    # Any columns not in `metaData` are moved to `featureData`
    #--------------------------------------#
    split_countdata = function(data) {
      data_tmp <- data.table::as.data.table(data)

      samples_in_metadata <- colnames(data_tmp) %in% private$.metaData[[ private$.sample_id ]]

      if (any(samples_in_metadata == FALSE)) {
        # Removes non-SAMPLE_IDs and adds them to the `featureData`
        private$.featureData <- data_tmp[, .SD, .SDcols = colnames(data_tmp)[!samples_in_metadata]]
        private$add_FEATURE_ID()
        cli::cli_alert_info("Placeholder {.field featureData} extracted from {.field countData}.")

        data_tmp[, (colnames(data_tmp)[!samples_in_metadata]) := NULL]
      }
      return(data_tmp)
    },
    # Checks & loads input table/filepath
    #--------------------------------------#
    check_table = function(data) {
      if (is.character(data) && length(data) == 1 && file.exists(data))
        return(data.table::fread(data, header = TRUE))

      if (inherits(data, "data.table") && !all(dim(data) == 0))
        return(data)

      if (is.data.frame(data) && !all(dim(data) == 0))
        return(data.table::as.data.table(data))

      return("Input must be an existing {.val filepath}, non-empty {.cls data.frame} or {.cls data.table}.")
    },
    # Checks & loads input matrix/filepath
    #--------------------------------------#
    check_matrix = function(data) {
      if (inherits(data, "sparseMatrix") && !all(dim(data) == 0))
        return(data)
        
      if ((is.matrix(data) || inherits(data, "denseMatrix")) && !all(dim(data) == 0))
        return(methods::as(data, "CsparseMatrix"))
        
      if (is.character(data) && length(data) == 1 && file.exists(data)) {
        dt <- data.table::fread(data, header = TRUE)

        # Removes rownames if present
        if (!is.null(dt$V1)) {
          dt_rownames <- dt$V1
          dt[, V1 := NULL]
        } else {
          dt_rownames <- NULL
        }
      } else if (inherits(data, "data.frame") || inherits(data, "data.table")) {
        dt <- data.table::as.data.table(data)
        dt_rownames <- NULL
      } else {
        return("Input must be an existing {.val filepath}, non-empty {.cls matrix}, {.cls Matrix}, {.cls data.frame} or {.cls data.table}.")
      }

      if (nrow(dt) == 0) {
        return("{.val filepath}, {.cls data.frame} or {.cls data.table} cannot be empty.")
      }

      # Only applies during data field initialisation.
      if (is.null(private$.countData) && is.null(private$.featureData)) {
        dt <- private$split_countdata(data = dt)
      }

      ## Replace strings to zero's
      cols <- names(dt)
      dt[, (cols) := lapply(.SD, function(x) {
            x <- as.numeric(as.character(x))
            x[is.na(x)] <- 0
            x
          }), .SDcols = cols
        ]

      if (inherits(dt, "data.table")) {
        # Return CsparseMatrix
        mat <- Matrix::Matrix(
          data = as.matrix(dt),
          dimnames = list(dt_rownames, colnames(dt))
        )
        return(methods::as(mat, "CsparseMatrix"))
      }
    }
  )
)
