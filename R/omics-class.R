#' Abstract omics class
#'
#' @description This is the abstract class 'omics', contains a variety of methods that are inherited and applied in the omics classes:
#' \link{metagenomics}, proteomics and metabolomics. 
#'
#' @details
#' Every class is created with the \link[R6]{R6Class} method. Methods are either public or private, and only the public components are inherited by other omic classes.
#' The omics class by default uses a \link[Matrix]{sparseMatrix} and \link[data.table]{data.table} data structures for quick and efficient data manipulation and returns the object by reference, same as the R6 class.
#' The method by reference is very efficient when dealing with big data.
#' @import R6 Matrix
#' @importFrom jsonlite toJSON
#' @importFrom jsonvalidate json_validate
#' @importFrom ape keep.tip
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
        cli::cli_abort("Data input must be {.cls data.table} like {.field metaData}.")
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
        cli::cli_abort("Data input must be {.cls data.table} like {.field featureData}.")
      }
    },
    #' @field countData A dense or sparse \link[Matrix]{Matrix}.
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
        cli::cli_abort("Data input must be {.cls Matrix} like {.field countData}.")
      }
    }
  ),
  public = list(
    #' @description
    #' Wrapper function that is inherited and adapted for each omics class.
    #' The omics classes requires a metadata samplesheet, that is validated by the metadata_schema.json.
    #' It requires a column `SAMPLE_ID` and optionally a `SAMPLEPAIR_ID` can be supplied. 
    #' The `SAMPLE_ID` will be used to link the metaData to the countData, and will act as the key during subsetting of other columns.
    #' To create a new object use [`new()`](#method-new) method. Do notice that the abstract class only checks if the metadata is valid!
    #' The `countData` and `featureData` will not be checked, these are handled by the sub-classes. 
    #' Using the omics class to load your data is not supported and still experimental.
    #' @param countData A path to an existing file or a dense/sparse \link[Matrix]{Matrix} format.
    #' @param featureData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @param metaData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @return A new `omics` object.
    #'
    initialize = function(countData = NULL, featureData = NULL, metaData = NULL) {
      #-------------------#
      ###   metaData    ###
      #-------------------#
      if (!is.null(metaData)) {
        duplicated_sample_ids <- FALSE
        private$.metaData <- private$check_table(metaData)
        self$validate()

        if (private$.valid_schema) {
          cli::cli_alert_success("{.field metaData} template passed the JSON validation.")

          #--------------------------------------------------------------------#
          ## Checking for duplicated sample identifiers
          #--------------------------------------------------------------------#

          cli::cli_alert_info("Checking for duplicated identifiers ..")
          duplicated_sample_idx <- duplicated(private$.metaData, by = private$.sample_id)
          duplicated_sample_ids <- any(duplicated_sample_idx)
          if (duplicated_sample_ids) {
            duplicated_sample_names <- unique(private$.metaData[[private$.sample_id]][duplicated_sample_idx])
            cli::cli_abort(
              "Found duplicated: {.val {duplicated_sample_names}}\
              \n Make sure {.arg SAMPLE_ID} column contains {.strong unique} identifiers!"
            )
          }
          #--------------------------------------------------------------------#
          ## Disable samplepair_id if not supplied
          #--------------------------------------------------------------------#
          if (!column_exists(private$.samplepair_id, private$.metaData))
            private$.samplepair_id <- NULL

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

        if (column_exists(private$.feature_id, private$.featureData)) {
          duplicated_feature_idx <- duplicated(private$.featureData, by = private$.feature_id)
          duplicated_feature_ids <- any(duplicated_feature_idx)

          if (duplicated_feature_ids) {
            duplicated_feature_names <- unique(private$.featureData[[private$.feature_id]][duplicated_feature_idx])
            cli::cli_abort(
              "Found duplicated: {.val {duplicated_feature_names}} \
              \n Make sure {.arg FEATURE_ID} column contains {.strong unique} identifiers!"
            )
          }

        } else {
          FEATURE_ID <- paste0("feature_", 1:nrow(private$.featureData))
          private$.featureData[, private$.feature_id := FEATURE_ID]
          data.table::setcolorder(
            x = private$.featureData,
            neworder = c(private$.feature_id, base::setdiff(colnames(private$.featureData), private$.feature_id))
          )
        }
        cli::cli_alert_success("{.field featureData} is loaded.")
      }

      #-------------------#
      ###   countData   ###
      #-------------------#
      if (!is.null(countData)) {
        private$.countData <- private$check_matrix(countData)
        cli::cli_alert_success("{.field countData} is loaded.")

        if (is.null(private$.featureData)) {
          private$add_featureData()
          cli::cli_alert_warning("Created placeholder {.field featureData}.")
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
    #' Validates an input metadata against the JSON schema. The metadata should look as follows and should not contain any empty spaces.
    #' For example; \code{'sample 1'} is not allowed, whereas \code{'sample1'} is allowed!
    #' 
    #' Acceptable column headers:
    #' * SAMPLE_ID (required)
    #' * SAMPLEPAIR_ID (optional)
    #' * CONTRAST_ (optional), used for [`autoFlow()`](#method-autoFlow).
    #' * VARIABLE_ (optional), not supported yet.
    #' 
    #' This function is used during the creation of a new object via [`new()`](#method-new) to validate the supplied metadata 
    #' via a filepath or existing \link[data.table]{data.table} or \link[base]{data.frame}.
    #' 
    #' @return None
    validate = function() {
      # Creates temporary json file from metadata
      tmp_json <- base::tempfile(fileext = ".json")

      json_data <- jsonlite::toJSON(
        private$.metaData,
        dataframe = "rows",
        pretty = TRUE,
        auto_unbox = TRUE
        )

      writeLines(json_data, tmp_json)

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
    #' @description
    #' Displays parameters of the omics class via stdout.
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
    #' @return object in place
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
    #' Upon creation of a new `omics` object a small backup of the original data is created.
    #' Since modification of the object is done by reference and duplicates are not made, it is possible to `reset` changes to the class.
    #' The methods from the abstract class \link{omics} also contains a private method to prevent any changes to the original object when using methods such as \code{ordination} \code{alpha_diversity} or \code{$DFE}.
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
    #' taxa$transform(log2)
    #'
    #' # resets
    #' taxa$reset()
    #'
    #' # An inbuilt reset function prevents unwanted modification to the taxa object.
    #' taxa$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
    #'
    #' @return object in place
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
    #' Remove NAs from `metaData` and updates the `countData`.
    #' @param column The column from where NAs should be removed, this can be either a wholenumbers or characters. Vectors are also supported.
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
    #' @return object in place
    removeNAs = function(column) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (all(is.wholenumber(column)) && length(column) <= length(colnames(private$.metaData)))
        column <- colnames(private$.metaData[column])

      if (!is.character(column))
        cli::cli_abort("{.val {column}} needs to be a character or an integer.")

      if (!column_exists(column, private$.metaData))
        cli::cli_abort("{.val {column}} does not exist in the {.field metaData} or one of the specified columns is completely empty!")

      ## MAIN
      #--------------------------------------------------------------------#
      private$.metaData <- na.omit(private$.metaData, cols = column)
      private$sync()
      invisible(self)
    },
    #' @description
    #' Feature subset (based on `featureData`), automatically applies data synchronization.
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in `featureData`.
    #' Only rows for which all conditions evaluate to TRUE are kept.
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
    #' @return object in place
    feature_subset = function(...) {
      # Replace all NAs by empty string
      features <- data.table::copy(private$.featureData)
      features[, names(features) := lapply(.SD, function(x) {
        if (is.character(x)) ifelse(is.na(x), "", x) else x
      })]

      rows_to_keep <- features[, ...]
      private$.featureData <- private$.featureData[rows_to_keep, ]
      private$.countData <- private$.countData[rows_to_keep, ]
      private$sync()
      self$print()
      invisible(self)
    },
    #' @description
    #' Sample subset (based on `metaData`), automatically applies synchronization.
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in `metaData`.
    #' Only rows for which all conditions evaluate to TRUE are kept.
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
    #' @return object in place
    sample_subset = function(...) {
      # set order of columns
      private$.countData <- private$.countData[, private$.metaData[[ private$.sample_id ]], drop = FALSE]
      # subset columns and rows
      rows_to_keep <- private$.metaData[, ...]
      private$.metaData <- private$.metaData[rows_to_keep, ]
      # NAs can occur in rows_to_keep, which then doesnt work on sparse Matrix.
      private$.countData <- private$.countData[, private$.metaData[[ private$.sample_id ]] ]
      private$sync()
      self$print()
      invisible(self)
    },
    #' @description
    #' Samplepair subset (based on `metaData`), automatically applies synchronization.
    #' @param num_unique_pairs An integer value to define the number of pairs to subset. The default is NULL, 
    #' meaning the maximum number of unique pairs will be used to subset the data. 
    #' Let's say you have three samples for each pair, then the `num_unique_pairs` will be set to 3.
    #' 
    #' @return object in place
    samplepair_subset = function(num_unique_pairs = NULL) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.null(num_unique_pairs) && !is.wholenumber(num_unique_pairs))
        cli::cli_abort("{.val {num_unique_pairs}} must contain integers!")

      ## MAIN
      #--------------------------------------------------------------------#

      counts <- private$.metaData[, .(unique_count = data.table::uniqueN(SAMPLE_ID)), by = SAMPLEPAIR_ID]

      if (is.null(num_unique_pairs)) {
        num_unique_pairs <- counts[, max(unique_count)]
      }

      private$.metaData <- private$.metaData[SAMPLEPAIR_ID %in% counts[unique_count == num_unique_pairs, SAMPLEPAIR_ID]]
      private$sync()
      self$print()
      invisible(self)
    },
    #' @description
    #' Agglomerates features by column, automatically applies synchronization.
    #' @param feature_rank A character value or vector of columns to aggregate from the `featureData`.
    #' @param feature_filter A character value or vector of characters to remove features via regex pattern.
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
    #' @return object in place
    feature_merge = function(feature_rank, feature_filter = NULL) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank))
        cli::cli_abort("{.val {feature_rank}} needs to be a character or vector containing characters")

      if (!is.null(feature_filter) && !is.character(feature_filter))
        cli::cli_abort("{.val {feature_filter}} needs to be a character or vector containing characters")

      if (!column_exists(feature_rank, private$.featureData))
        cli::cli_abort("{.val {feature_rank}} does not exist in {.field featureData}!")

      ## MAIN
      #--------------------------------------------------------------------#
      # creates a subset of unique feature rank, hashes combined for each unique rank
      counts <- data.table::data.table()
      counts[, (private$.feature_id) := rownames(private$.countData)]

      # Supports multiple features
      features <- data.table::copy(private$.featureData[private$.featureData[[ feature_rank[1] ]] != "", ])

      # set keys
      data.table::setkey(counts, FEATURE_ID)
      data.table::setkey(features, FEATURE_ID)

      # Create groups by ID
      grouped_ids <- features[, .(IDs = list(FEATURE_ID)), by = feature_rank]
      counts_glom <- Matrix::Matrix(0,
                                    nrow = nrow(grouped_ids),
                                    ncol = ncol(private$.countData),
                                    dimnames = list(NULL, colnames(private$.countData)),
                                    sparse = TRUE)

      # Populate sparse matrix by colsums of identical taxa
      for (i in 1:nrow(grouped_ids)) {
        ids <- grouped_ids$IDs[[i]]
        if (length(ids) == 1) {
          counts_glom[i, ] <- private$.countData[ids, ]
        } else {
          counts_glom[i, ] <- Matrix::colSums(private$.countData[ids, ])
        }
      }

      # Prepare final self-components
      private$.featureData <- base::unique(features, by = feature_rank)
      # Fetch first ID from each list
      grouped_ids$ID_first <- sapply(grouped_ids$IDs, `[[`, 1)
      # Reorder by matching IDs
      private$.featureData <- private$.featureData[ base::order(base::match(private$.featureData[[ private$.feature_id ]], grouped_ids$ID_first)) ]
      private$.countData <- counts_glom

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
    #' Performs transformation on the positive values from the `countData`.
    #' @param FUN A function such as \code{log2}, \code{log}
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
    #' obj$transform(log2)
    #'
    #' @return object in place
    transform = function(FUN) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!inherits(FUN, "function"))
        cli::cli_abort("{FUN} must be a function!")

      ## MAIN
      #--------------------------------------------------------------------#

      private$.countData@x <- FUN(private$.countData@x)
      invisible(self)
    },
    #' @description
    #' Relative abundance computation by column sums on the `countData`.
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
    #' obj$normalize()
    #'
    #' @return object in place
    normalize = function() {
      private$.countData@x <- private$.countData@x / rep(Matrix::colSums(private$.countData), base::diff(private$.countData@p))
      invisible(self)
    },
    #' @description
    #' Rank statistics based on `featureData`
    #' @details
    #' Counts the number of features identified for each column, for example in case of 16S metagenomics it would be the number of OTUs or ASVs on different taxonomy levels.
    #' @param feature_ranks A vector of characters or integers that match the `featureData`.
    #' @param unique A boolean value to display only unique entries in `feature_ranks`.
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
    #' plt <- obj$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
    #' plt
    #' @return A \link[ggplot2]{ggplot} object.
    #'
    rankstat = function(feature_ranks, unique = FALSE) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_ranks))
        cli::cli_abort("{.val {feature_ranks}} needs to be of character or integer type.")

      if (all(is.wholenumber(feature_ranks)) && length(feature_ranks) > length(colnames(private$.featureData)))
        feature_ranks <- colnames(private$.featureData[feature_ranks])

      if (!column_exists(feature_ranks, private$.featureData))
        cli::cli_abort("Specified {.val {feature_ranks}} do not exist in the {.field featureData}.")

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
      long_values <- data.table::melt(data = values,
                                      measure.vars = names(values),
                                      variable.name = "variable",
                                      value.name = "counts")

      # Sets order level of taxonomic ranks
      long_values[, variable := factor(variable, levels = base::rev(feature_ranks))]


      # Returns rankstat plot
      return(long_values %>%
               ggplot(mapping = aes(x = variable,
                                    y = counts)) +
               geom_col(fill = "grey",
                        colour = "grey15",
                        linewidth = 0.25) +
               coord_flip() +
               geom_text(mapping = aes(label = counts),
                         hjust = -0.1,
                         fontface = "bold") +
               ylim(0, max(long_values$counts)*1.10) +
               theme_bw() +
               labs(x = "Rank",
                    y = "Number of features classified"))
    },
    #' @description
    #' Alpha diversity based on \link{diversity}
    #' @param col_name A character variable from the `metaData`.
    #' @param metric An alpha diversity metric as input to \link{diversity}.
    #' @param group_by A column name to perform grouped statistical test in \link{diversity_plot} (default: NULL).
    #' @param Brewer.palID A character name for the palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link{colormap}.
    #' @param evenness A boolean wether to divide diversity by number of species, see \link[vegan]{specnumber}.
    #' @param paired A boolean value to perform paired analysis in \link[stats]{wilcox.test} and samplepair subsetting via [`samplepair_subset()`](#method-samplepair_subset)
    #' @param p.adjust.method A character variable to specify the p.adjust.method to be used, default is 'fdr'.
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
    #' plt <- obj$alpha_diversity(col_name = "treatment",
    #'                            metric = "shannon")
    #'
    #' @returns A list of components:
    #'  * `div` A \link[base]{data.frame} from \link{diversity}.
    #'  * `stats` A pairwise statistics from \link[rstatix]{pairwise_wilcox_test}.
    #'  * `plot` A \link[ggplot2]{ggplot} object.
    #' 
    #' @seealso \link{diversity_plot}
    alpha_diversity = function(col_name,
                               metric = c("shannon", "invsimpson", "simpson"),
                               Brewer.palID = "Set2",
                               group_by = NULL,
                               evenness = FALSE,
                               paired = FALSE,
                               p.adjust.method = "fdr") {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(col_name) && length(col_name) != 1) {
        cli::cli_abort("{.val {col_name}} must be a character and of length 1")
      } else if (!column_exists(col_name, private$.metaData)) {
        cli::cli_abort("The specified {.val {col_name}} does not exist in the {.field metaData}.")
      }

      if (!c(p.adjust.method %in% p.adjust.methods))
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

      # Remove NAs when col_name is specified
      if (!is.null(col_name))
        self$removeNAs(col_name)

      if (!is.null(group_by)) {
        combined_cols <- c(col_name, group_by)
      } else {
        combined_cols <- col_name
      }

      # Subset by samplepair completion
      if ( paired && !is.null(private$.samplepair_id) )
        self$samplepair_subset()

      # Alpha diversity based on 'metric'
      div <- data.table::data.table(diversity(x = private$.countData, metric=metric))
      div[, (combined_cols) := private$.metaData[, .SD, .SDcols = c(combined_cols)]]
      # Adjusts for evenness
      if (evenness) div$V1 <- div$V1 / log(vegan::specnumber(div$V1))

      # get colors
      colors <- colormap(private$.metaData, col_name, Brewer.palID)

      # Create and saves plots
      plot_list$data <- div
      diversity_plt <- diversity_plot(
        data = na.omit(div),
        values = "V1",
        col_name = col_name,
        group_by = group_by,
        palette = colors,
        method = metric,
        paired = paired,
        p.adjust.method = p.adjust.method
        )

      plot_list$stats <- as.data.frame(diversity_plt$stats)
      plot_list$plot <- diversity_plt$plot

      return(plot_list)
    },
    #' @description
    #' Creates a table most abundant compositional features. Also assigns a color blind friendly palette for visualizations.
    #' @param feature_rank A character variable in `featureData` to aggregate via [`feature_merge()`](#method-feature_merge).
    #' @param feature_filter A character or vector of characters to removes features by regex pattern.
    #' @param col_name Optional, a character or vector of characters to add to the final compositional data output.
    #' @param feature_top A wholenumber of the top features to visualize, the max is 15, due to a limit of palettes.
    #' @param normalize A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @param Brewer.palID A character name for the palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link{colormap}.
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
    #' result <- obj$composition(feature_rank = "Genus",
    #'                           feature_filter = c("uncultured"),
    #'                           feature_top = 10)
    #'
    #' plt <- composition_plot(data = result$data,
    #'                         palette = result$palette,
    #'                         feature_rank = "Genus")
    #'
    #' @returns A list of components:
    #'  * `data` A \link[data.table]{data.table} of feature compositions.
    #'  * `palette` A \link[stats]{setNames} palette from \link{colormap}.
    #' 
    #' @seealso \link{composition_plot}
    composition = function(feature_rank,
                           feature_filter = NULL,
                           col_name = NULL,
                           normalize = TRUE,
                           feature_top = c(10, 15),
                           Brewer.palID = "RdYlBu") {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.null(col_name)) {
        if (!is.character(col_name) && length(col_name) != 1) {
          cli::cli_abort("{.val {col_name}} must be a character and of length 1")
        } else if (!column_exists(col_name, private$.metaData)) {
          cli::cli_abort("The specified {.val {col_name}} does not exist in the {.field metaData}.")
        }
      }

      if (!is.wholenumber(feature_top)) {
        cli::cli_abort("{.val {feature_top}} must be an integer!")
      } else if (feature_top > 15) {
        cli::cli_alert_warning("The {.val {feature_top}} is set to an integer higher than 15.\n This may lead that colors are difficult to be distinguished.\n For color-blind people it is recommended to use a feature_top of maximum 15.")
      }

      if (!is.character(Brewer.palID) && length(Brewer.palID) != 1)
        cli::cli_abort("{.val {Brewer.palID}} must be a character and of length 1")

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

      # Normalizes sample counts
      if (normalize)
        self$normalize()

      # Agglomerate by feature_rank
      self$feature_merge(feature_rank = feature_rank, feature_filter = feature_filter)

      # Remove NAs when col_name is specified
      if (!is.null(col_name))
        self$removeNAs(col_name)

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
      if (!is.null(col_name)) {
        composition_final <- base::merge(final_long,
                                         private$.metaData[, .SD, .SDcols = c(private$.sample_id, col_name)],
                                         by = private$.sample_id,
                                         all = TRUE,
                                         allow.cartesian = TRUE) %>%
          unique()
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
    #' thus far supports 'bray', 'jaccard', 'cosine', 'manhattan', 'jsd' (jensen-shannon divergence), 'canberra' and 'unifrac' when a tree is provided via `treeData`, see [`distance()`](#method-distance).
    #' @param weighted A boolean value, to use abundances (\code{weighted = TRUE}) or absence/presence (\code{weighted=FALSE}) (default: TRUE).
    #' @param normalized A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @param threads A wholenumber, indicating the number of threads to use (Default: 1).
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
    #' @seealso \link{bray}, \link{canberra}, \link{cosine}, \link{jaccard}, \link{jsd}, \link{manhattan}, \link{unifrac}
    distance = function(metric, normalized = TRUE, weighted = TRUE, threads = 1) {

      ## Error handling
      #--------------------------------------------------------------------#
      OPTIONS <- c(
        "bray", "jaccard", "cosine", "manhattan",
        "jsd", "canberra", "unifrac"
        )

      if (!is.character(metric) && length(metric) != 1) {
        cli::cli_abort("{.val {metric}} needs to be a character with a length of 1")
      } else if (!metric %in% OPTIONS) {
        cli::cli_abort("{.val {metric}} is not a valid metric. \nValid options: {.val {OPTIONS}}")
      }

      if (!is.wholenumber(threads))
        cli::cli_abort("{.val {threads}} need to be an integer!")

      if (is.null(private$.treeData) && metric == "unifrac")
        cli::cli_abort("The specified {.val {metric}} is invalid since no {.field treeData} is supplied.")

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

      # Normalizes counts
      if (normalized)
        self$normalize()

      distmat <- switch(
        metric,
        "unifrac" = OmicFlow::unifrac(x = private$.countData, tree = private$.treeData, weighted=weighted, normalized=normalized, threads=threads),
        "manhattan" = OmicFlow::manhattan(x = private$.countData, weighted=weighted, threads=threads),
        "canberra" = OmicFlow::canberra(x = private$.countData, weighted=weighted, threads=threads),
        "jaccard" = OmicFlow::jaccard(x = private$.countData, weighted=weighted, threads=threads),
        "bray" = OmicFlow::bray(x = private$.countData, weighted=weighted, threads=threads),
        "jsd" = OmicFlow::jsd(x = private$.countData, weighted=weighted, threads=threads),
        "cosine" = OmicFlow::cosine(x = private$.countData, weighted=weighted, threads=threads)
      )

      return(distmat)
    },
    #' @description
    #' Ordination of `countData` with statistical testing.
    #' @param metric A dissimilarity or similarity metric to be applied on the `countData`, 
    #' thus far supports 'bray', 'jaccard', 'cosine', 'manhattan', 'jsd' (jensen-shannon divergence), 'canberra' and 'unifrac' when a tree is provided via `treeData`, see [`distance()`](#method-distance).
    #' @param method Ordination method, supports "pcoa" and "nmds", see \link[vegan]{wcmdscale}.
    #' @param distmat A custom distance matrix in either \link[stats]{dist} or \link[Matrix]{Matrix} format.
    #' @param group_by A character variable in `metaData` to be used for the \link{pairwise_adonis} or \link{pairwise_anosim} statistical test.
    #' @param weighted A boolean value, whether to compute weighted or unweighted dissimilarities (Default: TRUE).
    #' @param normalize A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @param threads A wholenumber, indicating the number of threads to use (Default: 1).
    #' @param perm_design A function that takes `metaData` and constructs a permutation design with \link[permute]{how} (default: NULL).
    #' @param perm A wholenumber, number of permutations to compare against the null hypothesis of \link[vegan]{adonis2} and \link[vegan]{anosim} (default: \code{perm=999}).
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
    #'                              weighted = TRUE,
    #'                              normalize = TRUE)
    #' pcoa_plots
    #'
    #' @returns A list of components:
    #'  * `distmat` A distance dissimilarity in \link[base]{matrix} format.
    #'  * `stats` A statistical test as a \link[base]{data.frame}.
    #'  * `pcs` principal components as a \link[base]{data.frame}.
    #'  * `scree_plot` A \link[ggplot2]{ggplot} object.
    #'  * `anova_plot` A \link[ggplot2]{ggplot} object.
    #'  * `scores_plot` A \link[ggplot2]{ggplot} object.
    #' 
    #' @seealso \link{ordination_plot}, \link{plot_pairwise_stats}, \link{pairwise_anosim}, \link{pairwise_adonis}
    ordination = function(metric = "bray",
                          method = c("pcoa", "nmds"),
                          group_by,
                          distmat = NULL,
                          weighted = TRUE,
                          normalize = TRUE,
                          threads = 1,
                          perm_design = NULL,
                          perm = 999) {

      ## Error handling
      #--------------------------------------------------------------------#
      if (!is.character(method) && length(method) != 1)
        cli::cli_abort("{.val {method}} needs to be a character with a length of 1")

      if (!is.character(group_by) && length(group_by) != 1) {
        cli::cli_abort("{.val {group_by}} needs to be a character with a length of 1")
      } else if (!column_exists(group_by, private$.metaData)) {
        cli::cli_abort("{.val {group_by}} does not exist in the metaData or is empty.")
      }

      if (!is.null(perm_design) && !is.function(perm_design))
        cli::cli_abort("perm_design must be a function.")

      if (!is.wholenumber(perm))
        cli::cli_abort("Permutations {.val {perm}} need to be an integer")

      if (!is.null(distmat) && (!inherits(distmat, "Matrix") && !inherits(distmat, "dist")))
        cli::cli_abort("{.arg distmat} need to be {.cls Matrix} or {.cls dist}")

      if (is.null(private$.treeData) && metric == "unifrac") {
        cli::cli_alert_warning("The specified {.val {metric}} is invalid since no tree is supplied.\n Switching to bray-curtis metric.")
        metric <- "bray"
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
          normalized = normalize,
          weighted = weighted,
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

      if (method == "pcoa") {
        # Normalisation of eigenvalues
        pcs$eig_norm <- unlist(lapply(pcs$eig, function(x) x / sum(pcs$eig) * 100))
        colnames(df_pcs_points) <- paste0("PC", 1:ncol(df_pcs_points))

      } else if (method == "nmds") {
        df_pcs_points[['stress']] <- pcs$stress
      }
      plot_list$pcs <- df_pcs_points

      # Adds relevant data
      df_pcs_points[, groups := private$.metaData[[ group_by ]] ]
      df_pcs_points[, samples := row.names(df_pcs_points) ]

      if (method == "pcoa") {
        # Scree plot of first 10 dimensions
        plot_list$scree_plot <- data.table::data.table(
          dims = seq(length(pcs$eig_norm[1:10])),
          dims.explained = pcs$eig_norm[1:10]
        ) %>%
          ggplot(mapping = aes(x = dims,
                               y = dims.explained)) +
          geom_col() +
          theme_bw() +
          scale_x_continuous(breaks=seq(1, 10, 1)) +
          scale_y_continuous(breaks=seq(0, 100, 10)) +
          labs(title = paste0("Screeplot of ", length(pcs$eig_norm[1:10])," PCs"),
               x = "Principal Components (PCs)",
               y = "dissimilarity explained [%]")

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
          col_name = "groups",
          pair=c("PC1", "PC2"),
          dist_explained = pcs$eig_norm,
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
          col_name = "groups",
          pair=c("MDS1", "MDS2"),
          dist_metric = metric
        )
      }

      return(plot_list)
    },
    #' @description
    #' Differential feature expression (DFE) using the \link{foldchange} for both paired and non-paired test.
    #' @param feature_rank A character or vector of characters in the `featureData` to aggregate via [`feature_merge()`](#method-feature_merge).
    #' @param feature_filter A character or vector of characters to remove features via regex pattern (Default: NULL).
    #' @param paired A boolean value, the paired is only applicable when a `SAMPLEPAIR_ID` column exists within the `metaData`. See \link[stats]{wilcox.test} and [`samplepair_subset()`](#method-samplepair_subset).
    #' @param condition.group A character variable of an existing column name in `metaData`, wherein the conditions A and B are located.
    #' @param condition_A A character value or vector of characters.
    #' @param condition_B A character value or vector of characters.
    #' @param pvalue.threshold A numeric value used as a p-value threshold to label and color significant features (Default: 0.05).
    #' @param logfold.threshold A numeric value used as a fold-change threshold to label and color significantly expressed features (Default: 0.06).
    #' @param abundance.threshold A numeric value used as an abundance threshold to size the scatter dots based on their mean abundance (default: 0.01).
    #' @param normalize A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
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
    #' unpaired <- obj$DFE(feature_rank = "Genus",
    #'                     paired = FALSE,
    #'                     condition.group = "treatment",
    #'                     condition_A = c("healthy"),
    #'                     condition_B = c("tumor"))
    #'
    #' @returns
    #'  * `dfe` A long \link[data.table]{data.table} table.
    #'  * `volcano_plot` A \link[ggplot2]{ggplot} object.
    #'
    #' @seealso \link{volcano_plot}, \link{foldchange}
    DFE = function(feature_rank,
                   feature_filter = NULL,
                   paired = FALSE,
                   normalize = TRUE,
                   condition.group,
                   condition_A,
                   condition_B,
                   pvalue.threshold = 0.05,
                   logfold.threshold = 0.06,
                   abundance.threshold = 0
                   ) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank) && length(feature_rank) != 1)
        cli::cli_abort("{.val {feature_rank}} needs to be a character with a length of 1")

      if (!is.character(condition.group) && length(condition.group) != 1) {
        cli::cli_abort("{.val {condition.group}} needs to be a character with a length of 1")
      } else if (!column_exists(condition.group, private$.metaData)) {
        cli::cli_abort("{.val {condition.group}} does not exist in the metaData or is empty.")
      }
      if (!is.character(condition_A))
        cli::cli_abort("{.val {condition_A}} needs to be a character.")

      if (!is.character(condition_B))
        cli::cli_abort("{.val {condition_B}} needs to be a character.")

      if (!is.numeric(pvalue.threshold))
        cli::cli_abort("{.val {pvalue.threshold}} need to be numeric.")

      if (!is.numeric(logfold.threshold))
        cli::cli_abort("{.val {logfold.threshold}} need to be numeric.")

      if (paired && is.null(private$.samplepair_id)) {
        cli::cli_alert_warning("Paired is set to {.val {paired}} but {.arg SAMPLEPAIR_ID} does not exist in the {.field metaData}.\n Differential feature analysis will continue now with paired set to {.val FALSE}!")
        paired <- FALSE
      }

      ## MAIN
      #--------------------------------------------------------------------#

      # Final output
      plot_list <- list()

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

      # normalization if applicable
      if (normalize)
        self$normalize()

      # Subset by missing values
      self$removeNAs(condition.group)

      # Subset by samplepair completion
      if (paired && !is.null(private$.samplepair_id))
        self$samplepair_subset()

      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_merge(feature_rank = feature_rank,
                        feature_filter = feature_filter)

      # Extract mean abundance
      abun <- as.matrix(Matrix::rowSums(private$.countData) / ncol(private$.countData))
      rownames(abun) <- private$.featureData[[ feature_rank ]]

      # Get data.table format abundances
      dt <- matrix_to_dtable(private$.countData)[, (feature_rank) := private$.featureData[[feature_rank]]]

      # Compute 2-fold expression based on (un)paired samples
      # Computes on equation of log2(A) - log2(B)
      # Supports multiple inputs for A and B.
      condition.labels <- data.table::setorderv(private$.metaData,
                                                cols = c(private$.sample_id, condition.group))[[ condition.group ]]

      # paired samples
      dfe <- foldchange(
        data = dt,
        condition_A = condition_A,
        condition_B = condition_B,
        paired = paired,
        condition_labels = condition.labels,
        feature_rank = feature_rank
        )

        #----------------------#
        # Visualization        #
        #----------------------#

        # Add abundance, and save data as output list
        dfe <- dfe[, "abun" := abun]
        plot_list$data <- dfe

        # Create & save volcano plot
        n_diff_columns <- sum(grepl("^Log2FC_", colnames(dfe)))

        plot_list$volcano_plot <- lapply(1:n_diff_columns, function(i) {
          volcano_plot(data = dfe,
                       logfold_col = paste0("Log2FC_", i),
                       pvalue_col = paste0("pvalue_", i),
                       feature_rank = feature_rank,
                       abundance_col = "abun",
                       pvalue.threshold = pvalue.threshold,
                       logfold.threshold = logfold.threshold,
                       abundance.threshold = abundance.threshold,
                       label_A = condition_A,
                       label_B = condition_B) +
            labs(
              subtitle = paste0(
                "Attribute: ", condition.group,
                ", test: ", ifelse(paired, "Wilcox signed rank test", "Mann-Whitney U test")
                )
            )
        })

      return(plot_list)
    },
    #' @description
    #' Automated Omics Analysis based on the `metaData`, see [`validate()`](#method-validate).
    #' For now only works with headers that start with prefix `CONTRAST_`. If the data is from the class `omics` or `proteomics` than FDR adjusted p-values are computed for the volcano plots.
    #' @param feature_contrast A character vector of feature columns in the `featureData` to aggregate via [`feature_merge()`](#method-feature_merge) (default: \code{"FEATURE_ID"}).
    #' @param feature_filter A character vector to filter unwanted features, (default: \code{NULL}).
    #' @param feature_ranks A character vector as input to [`rankstat()`](#method-rankstat) (default: \code{NULL}).
    #' @param distance_metrics A character vector specifying what (dis)similarity metrics to use (default: \code{c("unifrac")}).
    #' @param beta_div_table A path to an existing file or a dense/sparse \link[Matrix]{Matrix} format (default: \code{NULL}).
    #' @param alpha_div_table A path to pre-computed alpha diversity table, with columns: `alpha_div` (containing diversity values) and the same CONTRAST columns from `metaData` (default: \code{NULL}).
    #' @param normalize A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (default: \code{TRUE}).
    #' @param weighted A boolean value, whether to compute weighted or unweighted dissimilarities (default: \code{TRUE}).
    #' @param pvalue.threshold A numeric value, the p-value is used to include/exclude composition and foldchanges plots coming from alpha- and beta diversity analysis (default: 0.05).
    #' @param logfold.threshold A numeric value used as a fold-change threshold to label and color significantly expressed features, see [`DFE()`](#method-DFE) (Default: 1).
    #' @param abundance.threshold A numeric value used as an abundance threshold to size the scatter dots based on their mean abundance, see [`DFE()`](#method-DFE) (default: 0.01).
    #' @param perm A wholenumber, number of permutations to compare against the null hypothesis of \link[vegan]{adonis2} or \link[vegan]{anosim} (default: 999).
    #' @param threads Number of threads to use, only used in [`distance()`](#method-distance) when beta_div_table is not supplied (default: 1).
    #' @param report A boolean value to create a HTML markdown report (default: \code{FALSE}). If \code{FALSE} a nested list of the plots and data is returned.
    #' @param filename A character to name the HTML report to be saved in the current working directory (default: \code{paste0(getwd(), "/report.html")}). The \code{getwd()} is required for rmarkdown to save it in the right path.
    #' @importFrom patchwork plot_layout wrap_plots
    #' @return List of plots/data or rendered HTML report
    autoFlow = function(feature_contrast = "FEATURE_ID",
                        feature_filter = NULL,
                        feature_ranks = NULL,
                        distance_metrics = c("unifrac"),
                        beta_div_table = NULL,
                        alpha_div_table = NULL,
                        normalize = TRUE,
                        weighted = TRUE,
                        pvalue.threshold = 0.05,
                        logfold.threshold = 1,
                        abundance.threshold = 0.01,
                        perm = 999,
                        threads = 1,
                        report = TRUE,
                        filename = paste0(getwd(), "/report.html")
                      ) {
    ## Error handling
    #--------------------------------------------------------------------#

    if (!is.character(filename) && length(filename) != 1)
      cli::cli_abort("{.val {filename}} needs to be a character with a length of 1")
      
    if (!is.character(feature_contrast) && length(feature_contrast) != 1) {
      cli::cli_abort("{.val {feature_contrast}} needs to be a character with a length of 1")
    } else if (!column_exists(feature_contrast, private$.featureData)) {
      cli::cli_abort("{.val {feature_contrast}} does not exist in {.field featureData}!")
    }

    if (!is.null(beta_div_table) && !is.character(beta_div_table) && length(beta_div_table) != 1) {
      cli::cli_abort("{.arg beta_div_table} needs to be a character with a length of 1")
    
      if (!file.exists(beta_div_table))
        cli::cli_abort("{.arg beta_div_table} does not exists!")
    }

    if (!is.null(alpha_div_table) && !is.character(alpha_div_table) && length(alpha_div_table) != 1) {
      cli::cli_abort("{.arg alpha_div_table} needs to be a character with a length of 1")

      if (!file.exists(alpha_div_table))
        cli::cli_abort("{.arg alpha_div_table} does not exists!")
    }

    ## MAIN
    #--------------------------------------------------------------------#
    is_empty = function(obj) {
      if (length(obj) == 0) {
        return(NULL)
      } else {
        return(obj)
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
      if (!is.null(beta_div_table)) {
        beta_div_table <- private$check_matrix(filepath = beta_div_table)
        beta_div_table <- beta_div_table[private$.metaData[[private$.sample_id]], private$.metaData[[private$.sample_id]]]
      }

      # Load custom rarefraction alpha diversity table if supplied
      if (!is.null(alpha_div_table)) {
        alpha_div_table <- private$check_table(alpha_div_table)
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
        if (inherits(alpha_div_table, "data.table")) {
          res <- diversity_plot(
            data = alpha_div_table,
            values = "alpha_div",
            col_name = col_name,
            palette = colormap(dt_final, col_name, "Set2"),
            method = "custom"
            )
        } else {
          res <- tryCatch(
            {
              # Default attempt
              self$alpha_diversity(
                col_name = col_name,
                metric = "shannon",
                paired = ifelse(!is.null(private$.samplepair_id), TRUE, FALSE)
              )
            },
            error = function(e) {
              cli::cli_alert_warning("alpha_diversity with paired=TRUE failed. Retrying with paired=FALSE.")
              self$alpha_diversity(
                col_name = col_name,
                metric = "shannon",
                paired = FALSE
              )
            }
          )
        }
        
        ## Save plots & data
        alpha_div_plots[[i]] <- res$plot
        alpha_div_data[[i]] <- list(data = res$data, stats = res$stats)

        ### Identify significant groups for composition plots & volcano plots
        signif_pairs <- res$stats[res$stats$p.adj < pvalue.threshold, ][c("group1", "group2")]
        if (nrow(signif_pairs) > 0)
          conditions <- signif_pairs
          
        #--------------------------------------------------------------------#
        ## Beta diversity
        #--------------------------------------------------------------------#
        
        for (j in 1:metrics_nrow) {
          if (inherits(beta_div_table, "Matrix")) {
            res <- self$ordination(
              distmat = beta_div_table,
              method = "pcoa",
              perm = perm,
              group_by = col_name
              )
          } else {
            res <- self$ordination(
              metric = distance_metrics[j],
              method = "pcoa",
              group_by = col_name,
              normalize = normalize,
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
          if (inherits(beta_div_table, "Matrix")) {
            res <- self$ordination(
              distmat = beta_div_table,
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
              normalize = normalize,
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
          # Creates composition long table
          res <- self$composition(
            feature_rank = feature_contrast[j],
            feature_filter = feature_filter,
            feature_top = 15,
            normalize = normalize,
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
          
          if (!is.null(conditions) && nrow(conditions) > 0) {

            dfe <- tryCatch(
              {
              # Default attempt
              self$DFE(
                feature_rank = feature_contrast[j],
                feature_filter = feature_filter,
                paired = ifelse(!is.null(private$.samplepair_id), TRUE, FALSE),
                normalize = normalize,
                condition.group = col_name,
                condition_A = c(conditions$group1),
                condition_B = c(conditions$group2),
                pvalue.threshold = pvalue.threshold,
                abundance.threshold = abundance.threshold,
                logfold.threshold = logfold.threshold
                )
              },
              error = function(e) {
                cli::cli_alert_warning("DFE with paired=TRUE failed. Retrying with paired=FALSE.")
                self$DFE(
                  feature_rank = feature_contrast[j],
                  feature_filter = feature_filter,
                  paired = FALSE,
                  normalize = normalize,
                  condition.group = col_name,
                  condition_A = c(conditions$group1),
                  condition_B = c(conditions$group2),
                  pvalue.threshold = pvalue.threshold,
                  abundance.threshold = abundance.threshold,
                  logfold.threshold = logfold.threshold
                  )
              }
            )
            if (class(self)[1] %in% c("omics", "proteomics")) {
              dfe$data$p.adj <- p.adjust(p = dfe$data$pvalue_1, method = "fdr")
              dfe$volcano_plot <- volcano_plot(
                data = dfe$data,
                logfold_col = "Log2FC_1",
                pvalue_col = "p.adj",
                feature_rank = feature_contrast[j],
                abundance_col = "abun",
                label_A = conditions$group1,
                label_B = conditions$group2,
                pvalue.threshold = pvalue.threshold,
                abundance.threshold = abundance.threshold,
                logfold.threshold = logfold.threshold
              )
            }
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
          common_samples <- base::intersect(private$.metaData[[ private$.sample_id ]], colnames(private$.countData))
          private$.countData <- private$.countData[, common_samples, drop = FALSE]
          private$.metaData <- private$.metaData[private$.metaData[[ private$.sample_id ]] %in% common_samples, ]
        }
      }

      if (!is.null(private$.featureData)) {
        if (!column_exists(private$.feature_id, private$.featureData))
          cli::cli_abort("{private$.feature_id} doesn't exist in {.field featureData}.")

        private$.featureData <- private$check_table(private$.featureData)
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
        cli::cli_alert_warning("Placeholder {.field featureData} created.")
      }
    },
    removeZeros = function() {
      keep_cols <- Matrix::colSums(private$.countData) > 0
      keep_rows <- Matrix::rowSums(private$.countData) > 0

      private$.countData <- private$.countData[keep_rows, keep_cols]
      private$.metaData <- private$.metaData[keep_cols, ]
      private$.featureData <- private$.featureData[keep_rows]

      if (!is.null(private$.treeData))
        private$.treeData <- ape::keep.tip(private$.treeData, private$.featureData[[ private$.feature_id ]])
    },
    add_featureData = function() {
      private$.featureData <- data.table::data.table()
      countData_with_rownames <- rownames(private$.countData)

      if (is.null(countData_with_rownames)) {
        FEATURE_ID <- paste0("feature_", 1:nrow(private$.countData))
        private$.featureData <- private$.featureData[, (private$.feature_id) := FEATURE_ID]
        rownames(private$.countData) <- FEATURE_ID
      } else {
        private$.featureData <- private$.featureData[, (private$.feature_id) := countData_with_rownames]
      }          
    },
    # Checks & loads input table/filepath
    #--------------------------------------#
    check_table = function(data) {
    if (is.character(data) && length(data) == 1 && file.exists(data))
      return(data.table::fread(data, header = TRUE))

    if (inherits(data, "data.table"))
      return(data)

    if (is.data.frame(data))
      return(data.table::as.data.table(data))

    cli::cli_abort("Input must be an existing {.val filepath}, {.cls data.frame} or {.cls data.table}.")
  },

  # Checks & loads input matrix/filepath
  #--------------------------------------#
  check_matrix = function(data) {
    if (is.character(data) && length(data) == 1 && file.exists(data)) {
      dt <- data.table::fread(data, header = TRUE)
      # Change character values to numeric
      for (col in names(dt)) {
        dt[is.na(get(col)), (col) := 0]
        dt[get(col) == "", (col) := 0]
      }

      # Removes rownames if present
      if (!is.null(dt$V1)) {
        dt_rownames <- dt$V1
        dt[, V1 := NULL]
      } else {
        dt_rownames <- NULL
      }
      # Convert to matrix format
      mat <- Matrix::Matrix(
        data = as.matrix(dt),
        dimnames = list(dt_rownames, colnames(dt))
      )
      
      # Return CsparseMatrix
      return(as(mat, "CsparseMatrix"))
    }

    if (inherits(data, "sparseMatrix"))
      return(data)

    if (is.matrix(data) || inherits(data, "denseMatrix"))
      return(as(data, "CsparseMatrix"))
      
    cli::cli_abort("Input must be an existing {.val filepath}, {.cls matrix} or {.cls Matrix}.")
    }
  )
)
