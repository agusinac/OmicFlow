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
#' @export

omics <- R6::R6Class(
  classname = "omics",
  cloneable = FALSE,
  public = list(
    #' @field countData A path to an existing file, data.table or data.frame.
    countData = NULL,

    #' @field featureData A path to an existing file, data.table or data.frame.
    featureData = NULL,

    #' @field metaData A path to an existing file, data.table or data.frame.
    metaData = NULL,

    #' @field .valid_schema Boolean value for schema validation via JSON
    .valid_schema = NULL,

    #' @field .feature_id A character, default name for the feature identifiers.
    .feature_id = "FEATURE_ID",

    #' @field .sample_id A character, default name for the sample identifiers.
    .sample_id = "SAMPLE_ID",

    #' @field .samplepair_id A character, default name for the sample pair identifiers.
    .samplepair_id = "SAMPLEPAIR_ID",

    #' @description
    #' Wrapper function that is inherited and adapted for each omics class.
    #' The omics classes requires a metadata samplesheet, that is validated by the metadata_schema.json.
    #' It requires a column `SAMPLE_ID` and optionally a `SAMPLEPAIR_ID` or `FEATURE_ID` can be supplied. 
    #' The `SAMPLE_ID` will be used to link the metaData to the countData, and will act as the key during subsetting of other columns.
    #' To create a new object use [`new()`](#method-new) method. Do notice that the abstract class only checks if the metadata is valid!
    #' The `countData` and `featureData` will not be checked, these are handles by the sub-classes. 
    #' Using the omics class to load your data is not supported and still experimental.
    #' @param countData A path to an existing file, data.table or data.frame.
    #' @param featureData A path to an existing file, data.table or data.frame.
    #' @param metaData A path to an existing file, data.table or data.frame.
    #' @return A new `omics` object.
    #'
    initialize = function(countData = NULL, featureData = NULL, metaData = NULL) {
      #-------------------#
      ###   metaData    ###
      #-------------------#
      if (!is.null(metaData)) {
        self$metaData <- private$check_table(metaData)
        self$validate()

        if (self$.valid_schema) {
          cli::cli_alert_success("Metadata template passed the JSON validation.")

          self$metaData <- self$metaData[, lapply(.SD, function(x) ifelse(x == "", NA, x)),
                                         .SDcols = colnames(self$metaData)]

          colnames(self$metaData) <- gsub("\\s+", "_", colnames(self$metaData))

          #--------------------------------------------------------------------#
          ## Checking for duplicated sample and feature identifiers
          #--------------------------------------------------------------------#

          cli::cli_alert_info("Checking for duplicated identifiers ..")

          duplicated_sample_ids <- any(duplicated(self$metaData, by = self$.sample_id))

          if (column_exists(self$.feature_id, self$metaData)) {
            duplicated_feature_ids <- any(duplicated(self$metaData, by = self$.feature_id))
          } else {
            duplicated_feature_ids <- FALSE
          }

          if (duplicated_sample_ids) {
            cli::cli_abort("Found duplicated SAMPLE_ID, make sure SAMPLE_ID column contains unique identifiers!")
          } else if (duplicated_feature_ids) {
            cli::cli_abort("Found duplicated FEATURE_ID, make sure FEATURE_ID column contains unique identifiers!")
          }

          #--------------------------------------------------------------------#
          ## Disable samplepair_id if not supplied
          #--------------------------------------------------------------------#
          if (!column_exists(self$.samplepair_id, self$metaData))
            self$.samplepair_id <- NULL

        } else {
          errors <- attr(self$.valid_schema, "errors")
          cli::cli_abort(
            "JSON validation failed: \n{ paste(errors$message, collapse = '\n')}"
            )
        }

      } else {
        cli::cli_abort(
          "metaData cannot be empty, please provide a data.frame, data.table or filepath"
        )
      }

      #-------------------#
      ###  featureData  ###
      #-------------------#
      if (!is.null(featureData)) {
        self$featureData <- private$check_table(featureData)

        if (column_exists(self$.feature_id, self$metaData)) {
          FEATURE_ID <- self$metaData[[self$.feature_id]]
        } else {
          FEATURE_ID <- paste0("feature_", rownames(self$featureData))
        }

        self$featureData[, (self$.feature_id) := FEATURE_ID]
        colnames(self$featureData) <- gsub("\\s+", "_", colnames(self$featureData))

        cli::cli_alert_success("featureData is loaded.")
      }

      #-------------------#
      ###   countData   ###
      #-------------------#
      if (!is.null(countData)) {
        self$countData <- private$check_matrix(countData)
        rownames(self$countData) <- self$featureData$FEATURE_ID
        cli::cli_alert_success("countData is loaded.")
      }

      # There should be an interal metadata template check to make sure all headers are correct.
      # Should also include to check for missing data and alert the user!
      # Current example to be used in the future
      base::tryCatch(
        { self$countData <- self$countData[, self$metaData[[ self$.sample_id ]], drop = FALSE] },
        error = function(e) {
          cli::cli_abort(c(
            "Error occured during countData subsetting by metaData:",
            "i" = "The error message was: {e$message}"
          ))
        }
      )

    },
    #' @description
    #' Validates an input metadata against the JSON schema. The metadata should look as follows and should not contain any empty spaces.
    #' For example; \code{'sample 1'} is not allowed, whereas \code{'sample1'} is allowed!
    #' 
    #' Acceptable column headers:
    #' * SAMPLE_ID (required)
    #' * SAMPLEPAIR_ID (optional)
    #' * FEATURE_ID (optional)
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
        self$metaData,
        dataframe = "rows",
        pretty = TRUE,
        auto_unbox = TRUE
        )

      writeLines(json_data, tmp_json)

      self$.valid_schema <- jsonvalidate::json_validate(
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
    #' Removes empty (zero) values by row and column from the `countData`.
    #' This method is performed automatically during subsetting of the object.
    #' @examples
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
    #' 
    #' obj$removeZeros()
    #' 
    #' @return object in place
    removeZeros = function() {
      # Remove empty samples (columns)
      keep_cols <- Matrix::colSums(self$countData) > 0

      # Remove empty species (rows)
      keep_rows <- Matrix::rowSums(self$countData) > 0

      # Creates new countData instance
      self$countData <- self$countData[keep_rows, keep_cols]
      self$featureData <- self$featureData[keep_rows]
      invisible(self)
    },
    #' @description
    #' Remove NAs from `metaData` and updates the `countData`.
    #' @param column The column from where NAs should be removed, this can be either a wholenumbers or characters. Vectors are also supported.
    #' @examples
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
    #' 
    #' obj$removeNAs(column = "treatment")
    #' 
    #' @return object in place
    removeNAs = function(column) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (all(is.wholenumber(column)) && length(column) <= length(colnames(self$metaData)))
        column <- colnames(self$metaData[column])

      if (!is.character(column))
        cli::cli_abort("{column} needs to be a character or an integer.")

      if (!column_exists(column, self$metaData))
        cli::cli_abort("{column} do not exist in the metaData or one of the specified columns is completely empty!")

      ## MAIN
      #--------------------------------------------------------------------#

      self$metaData <- na.omit(self$metaData, cols = column)
      self$countData <- self$countData[, self$metaData[[ self$.sample_id ]]]
      invisible(self)
    },
    #' @description
    #' Feature subset (based on `featureData`), automatically applies [`removeZeros()`](#method-removeZeros).
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in `featureData`.
    #' Only rows for which all conditions evaluate to TRUE are kept.
    #' @examples
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
    #' 
    #' obj$feature_subset(Genus == "Streptococcus")
    #' 
    #' @return object in place
    feature_subset = function(...) {
      # Replace all NAs by empty string
      features <- data.table::copy(self$featureData)
      features[, names(features) := lapply(.SD, function(x) {
        if (is.character(x)) ifelse(is.na(x), "", x) else x
      })]

      rows_to_keep <- features[, ...]
      self$featureData <- self$featureData[rows_to_keep, ]
      self$countData <- self$countData[rows_to_keep, ]
      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Sample subset (based on `metaData`), automatically applies [`removeZeros()`](#method-removeZeros).
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in `metaData`.
    #' Only rows for which all conditions evaluate to TRUE are kept.
    #' @examples
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
    #' 
    #' obj$sample_subset(treatment == "tumor")
    #'
    #' @return object in place
    sample_subset = function(...) {
      # set order of columns
      self$countData <- self$countData[, self$metaData[[ self$.sample_id ]], drop = FALSE]
      # subset columns and rows
      rows_to_keep <- self$metaData[, ...]
      self$metaData <- self$metaData[rows_to_keep, ]
      # NAs can occur in rows_to_keep, which then doesnt work on sparse Matrix.
      self$countData <- self$countData[, self$metaData[[ self$.sample_id ]] ]
      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Samplepair subset (based on `metaData`), automatically applies [`removeZeros()`](#method-removeZeros).
    #' @param num_unique_pairs An integer value to define the number of pairs to subset. The default is NULL, 
    #' meaning the maximum number of unique pairs will be used to subset the data. 
    #' Let's say you have three samples for each pair, then the `num_unique_pairs` will be set to 3.
    #' 
    #' @return object in place
    samplepair_subset = function(num_unique_pairs = NULL) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.null(num_unique_pairs) && !is.wholenumber(num_unique_pairs))
        cli::cli_abort("{num_unique_pairs} must contain integers!")

      ## MAIN
      #--------------------------------------------------------------------#

      counts <- self$metaData[, .(unique_count = data.table::uniqueN(SAMPLE_ID)), by = SAMPLEPAIR_ID]

      if (is.null(num_unique_pairs)) {
        num_unique_pairs <- counts[, max(unique_count)]
      }

      self$metaData <- self$metaData[SAMPLEPAIR_ID %in% counts[unique_count == num_unique_pairs, SAMPLEPAIR_ID]]
      self$countData <- self$countData[, self$metaData[[self$.sample_id]] ]
      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Agglomerates features by column, automatically applies [`removeZeros()`](#method-removeZeros).
    #' @param feature_rank A character value or vector of columns to aggregate from the `featureData`.
    #' @param feature_filter A character value or vector of characters to remove features via regex pattern.
    #' @examples
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
    #' 
    #' obj$feature_merge(feature_rank = c("Kingdom", "Phylum"))
    #' obj$feature_merge(feature_rank = "Genus", feature_filter = c("uncultured", "metagenome"))
    #'
    #' @return object in place
    feature_merge = function(feature_rank, feature_filter = NULL) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank))
        cli::cli_abort("{feature_rank} needs to be a character or vector containing characters")

      if (!is.null(feature_filter) && !is.character(feature_filter))
        cli::cli_abort("{feature_filter} needs to be a character or vector containing characters")

      if (!column_exists(feature_rank, self$featureData))
        cli::cli_abort("{feature_rank} does not exist in featureData!")

      ## MAIN
      #--------------------------------------------------------------------#

      # creates a subset of unique feature rank, hashes combined for each unique rank
      counts <- data.table::data.table("FEATURE_ID" = rownames(self$countData))

      # Supports multiple features
      features <- data.table::copy(self$featureData[self$featureData[[ feature_rank[1] ]] != "", ])

      # set keys
      data.table::setkey(counts, FEATURE_ID)
      data.table::setkey(features, FEATURE_ID)

      # Create groups by ID
      grouped_ids <- features[, .(IDs = list(FEATURE_ID)), by = feature_rank]
      counts_glom <- Matrix::Matrix(0,
                                    nrow = nrow(grouped_ids),
                                    ncol = ncol(self$countData),
                                    dimnames = list(NULL, colnames(self$countData)),
                                    sparse = TRUE)

      # Populate sparse matrix by colsums of identical taxa
      for (i in 1:nrow(grouped_ids)) {
        ids <- grouped_ids$IDs[[i]]
        if (length(ids) == 1) {
          counts_glom[i, ] <- self$countData[ids, ]
        } else {
          counts_glom[i, ] <- Matrix::colSums(self$countData[ids, ])
        }
      }

      # Prepare final self-components
      self$featureData <- base::unique(features, by = feature_rank)
      # Fetch first ID from each list
      grouped_ids$ID_first <- sapply(grouped_ids$IDs, `[[`, 1)
      # Reorder by matching IDs
      self$featureData <- self$featureData[ base::order(base::match(self$featureData$FEATURE_ID, grouped_ids$ID_first)) ]
      self$countData <- counts_glom

      # Replaces strings matching feature_filter with NAs
      if (!is.null(feature_filter)) {
        regex_pattern <- paste(feature_filter, collapse = "|")
        for (col in feature_rank) {
          self$featureData[
            grepl(regex_pattern, get(col), ignore.case = TRUE),
            (col) := NA_character_
          ]
        }
      }

      # Clean up featureData
      empty_strings <- !is.na(self$featureData[[ feature_rank[1] ]])
      self$featureData <- self$featureData[empty_strings, ]
      self$countData <- self$countData[empty_strings, ]
      rownames(self$countData) <- self$featureData$FEATURE_ID

      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Performs transformation on the positive values from the `countData`.
    #' @param FUN A function such as \code{log2}, \code{log}
    #' @examples
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
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

      self$countData@x <- FUN(self$countData@x)
      invisible(self)
    },
    #' @description
    #' Relative abundance computation by column sums on the `countData`.
    #' @examples
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
    #' 
    #' obj$normalize()
    #'
    #' @return object in place
    normalize = function() {
      self$countData@x <- self$countData@x / rep(Matrix::colSums(self$countData), base::diff(self$countData@p))
      invisible(self)
    },
    #' @description
    #' Rank statistics based on `featureData`
    #' @details
    #' Counts the number of features identified for each column, for example in case of 16S metagenomics it would be the number of OTUs or ASVs on different taxonomy levels.
    #' @param feature_ranks A vector of characters or integers that match the `featureData`.
    #' @examples
    #' library(ggplot2)
    #' 
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
    #' 
    #' plt <- obj$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"))
    #' plt
    #' @return A \link[ggplot2]{ggplot} object.
    #'
    rankstat = function(feature_ranks) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_ranks))
        cli::cli_abort("{feature_ranks} needs to be of character or integer type.")

      if (all(is.wholenumber(feature_ranks)) && length(feature_ranks) > length(colnames(self$featureData)))
        feature_ranks <- colnames(self$featureData[feature_ranks])

      if (!column_exists(feature_ranks, self$featureData))
        cli::cli_abort("Specified {feature_ranks} do not exist in the featureData.")

      ## MAIN
      #--------------------------------------------------------------------#

      # Counts number of ASVs without empty values
      values <- self$featureData[, lapply(.SD, function(x) sum(!is.na(x) & x != "")), .SDcols = !c(self$.feature_id)][, .SD, .SDcols = feature_ranks]

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
    #' @param Brewer.palID A character name for the palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link{colormap}.
    #' @param evenness A boolean wether to divide diversity by number of species, see \link[vegan]{specnumber}.
    #' @param paired A boolean value to perform paired analysis in \link[stats]{wilcox.test} and samplepair subsetting via [`samplepair_subset()`](#method-samplepair_subset)
    #' @param p.adjust.method A character variable to specify the p.adjust.method to be used, default is 'fdr'.
    #' @examples
    #' library(ggplot2)
    #' 
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
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
                               evenness = FALSE,
                               paired = FALSE,
                               p.adjust.method = "fdr") {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(col_name) && length(col_name) != 1) {
        cli::cli_abort("{col_name} must be a character and of length 1")
      } else if (!column_exists(col_name, self$metaData)) {
        cli::cli_abort("The specified {col_name} does not exist in the metaData.")
      }

      if (!c(p.adjust.method %in% p.adjust.methods))
        cli::cli_abort("Specified {p.adjust.method} is not valid. \nValid options: {p.adjust.methods}")

      ## MAIN
      #--------------------------------------------------------------------#

      # OUTPUT: Plot list
      plot_list <- list()

      # Save omics class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )
      
      # Restores omics class components
      on.exit(private$tmp_restore(), add = TRUE)

      # Remove NAs when col_name is specified
      if (!is.null(col_name))
        self$removeNAs(col_name)

      # Subset by samplepair completion
      if ( paired && !is.null(self$.samplepair_id) )
        self$samplepair_subset()

      # Alpha diversity based on 'metric'
      div <- data.table::data.table(diversity(x = self$countData, metric=metric))
      div[, (col_name) := self$metaData[, .SD, .SDcols = c(col_name)]]
      # Adjusts for evenness
      if (evenness) div$V1 <- div$V1 / log(vegan::specnumber(div$V1))

      # get colors
      colors <- colormap(self$metaData, col_name, Brewer.palID)

      # Create and saves plots
      plot_list$data <- div
      diversity_plt <- diversity_plot(
        data = na.omit(div),
        values = "V1",
        col_name = col_name,
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
    #' @importFrom viridis viridis
    #' @examples
    #' library(ggplot2)
    #' 
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
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
          cli::cli_abort("{col_name} must be a character and of length 1")
        } else if (!column_exists(col_name, self$metaData)) {
          cli::cli_abort("The specified {col_name} does not exist in the metaData.")
        }
      }

      if (!is.wholenumber(feature_top)) {
        cli::cli_abort("{feature_top} must be an integer!")
      } else if (feature_top > 15) {
        cli::cli_alert_warning("The {feature_top} is set to an integer higher than 15.\n This may lead that colors are difficult to be distinguished.\n For color-blind people it is recommended to use a feature_top of maximum 15.")
      }

      if (!is.character(Brewer.palID) && length(Brewer.palID) != 1)
        cli::cli_abort("{Brewer.palID} must be a character and of length 1")

      ## MAIN
      #--------------------------------------------------------------------#

      # Copies object to prevent modification of omics class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )
      
      # Restores omics class components
      on.exit(private$tmp_restore(), add = TRUE)

      # Normalizes sample counts
      if (normalize)
        self$normalize()

      # Agglomerate by feature_rank
      self$feature_merge(feature_rank = feature_rank, feature_filter = feature_filter)

      # Remove NAs when col_name is specified
      if (!is.null(col_name))
        self$removeNAs(col_name)

      # Convert sparse matrix to data.table
      counts <- sparse_to_dtable(self$countData)

      # Fetch unfiltered and filtered features
      dt <- counts[, (feature_rank) := self$featureData[[feature_rank]]]

      # Create row_sums
      dt[, row_sum := rowSums(.SD), .SDcols = !c(feature_rank)]

      # Orders by row_sum in descending order
      data.table::setorder(dt, -row_sum)

      # Subset taxa for visualization
      final_dt <- rbind(dt[1:feature_top][, .SD, .SDcols = !c("row_sum")],
                        dt[(feature_top+1):nrow(dt)][, lapply(.SD, function(x) sum(x)),
                                                                 .SDcols = !c(feature_rank, "row_sum")],
                        fill = TRUE)
      final_dt[nrow(final_dt), (feature_rank)] <- "Other"

      # Creates palette
      df_taxa_len <- length(final_dt[[feature_rank]])
      if (Brewer.palID == FALSE) {
        chosen_palette <- viridis(df_taxa_len - 1)
      } else if (df_taxa_len-1 <= 15 && df_taxa_len-1 > 10) {
        chosen_palette <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                            "#920000","#924900","#db6d00","#24ff24","#ffff6d")[1:df_taxa_len-1]
      } else {
        chosen_palette <- RColorBrewer::brewer.pal(df_taxa_len-1, Brewer.palID)
      }
      taxa_colors_ordered <- stats::setNames(c(chosen_palette, "lightgrey"), final_dt[[feature_rank]])

      # Pivoting in long table and factoring feature ranke
      final_long <- data.table::melt(final_dt,
                                     id.vars = c(feature_rank),
                                     variable.factor = FALSE,
                                     value.factor = TRUE)
      # Rename colnames for merge step
      colnames(final_long) <- c(feature_rank, self$.sample_id, "value")

      # Adds metadata columns by user input
      if (!is.null(col_name)) {
        composition_final <- base::merge(final_long,
                                         self$metaData[, .SD, .SDcols = c(self$.sample_id, col_name)],
                                         by = self$.sample_id,
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
    #' Ordination of `countData` with statistical testing.
    #' @param metric A dissimilarity or similarity metric to be applied on the `countData`, 
    #' thus far supports 'bray', 'jaccard' and 'unifrac' when a tree is provided via `treeData`, see \link[rbiom]{bdiv_distmat}.
    #' @param method Ordination method, supports "pcoa" and "nmds", see \link[vegan]{wcmdscale}.
    #' @param distmat A custom distance matrix in either \link[stats]{dist} or \link[Matrix]{Matrix} format.
    #' @param group_by A character variable in `metaData` to be used for the \link{pairwise_adonis} or \link{pairwise_anosim} statistical test.
    #' @param weighted A boolean value, whether to compute weighted or unweighted dissimilarities (Default: TRUE).
    #' @param normalize A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @param cpus A wholenumber, indicating the number of processes to spawn (Default: 1) in \link[rbiom]{bdiv_distmat}.
    #' @param perm A wholenumber, number of permutations to compare against the null hypothesis of \link[vegan]{adonis2} and \link[vegan]{anosim} (default: \code{perm=999}).
    #' @importFrom purrr map
    #' @importFrom rbiom bdiv_distmat
    #' @importFrom slam as.simple_triplet_matrix
    #' @examples
    #' library(ggplot2)
    #' 
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
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
    ordination = function(metric = c("bray", "jaccard", "unifrac"),
                          method = c("pcoa", "nmds"),
                          group_by,
                          distmat = NULL,
                          weighted = TRUE,
                          normalize = TRUE,
                          cpus = 1,
                          perm = 999) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(metric) && length(metric) != 1)
        cli::cli_abort("{metric} needs to be a character with a length of 1")

      if (!is.character(method) && length(method) != 1)
        cli::cli_abort("{method} needs to be a character with a length of 1")

      if (!is.character(group_by) && length(group_by) != 1) {
        cli::cli_abort("{group_by} needs to be a character with a length of 1")
      } else if (!column_exists(group_by, self$metaData)) {
        cli::cli_abort("{group_by} does not exist in the metaData or is empty.")
      }

      if (!is.wholenumber(cpus))
        cli::cli_abort("{cpus} need to be an integer!")

      if (!is.wholenumber(perm))
        cli::cli_abort("Permutations {perm} need to be an integer")

      if (!is.null(distmat) && (!inherits(distmat, "Matrix") && !inherits(distmat, "dist")))
        cli::cli_abort("custom distance matrix (distmat) need to be of class Matrix or dist")

      if (is.null(self$treeData) && metric == "unifrac") {
        cli::cli_alert_warning("The specified {metric} is invalid since no tree is supplied.\n Switching to bray-curtis metric.")
        metric <- "bray"
      }

      ## MAIN
      #--------------------------------------------------------------------#

      # Copies object to prevent modification of omics class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )
      
      # Restores omics class components
      on.exit(private$tmp_restore(), add = TRUE)

      # Subset by missing values
      self$removeNAs(group_by)
      if (inherits(distmat, "Matrix")) {
        distmat <- distmat[self$metaData[[ self$.sample_id ]], self$metaData[[ self$.sample_id ]]]
        distmat <- as.dist(distmat)
      }

      # Creates a list of plots
      plot_list <- list()

      # Normalizes counts
      if (normalize)
        self$normalize()

      # Requires rownames to contain same labels as tree
      if (is.null(distmat)) {
        counts <- slam::as.simple_triplet_matrix(self$countData)
        rownames(counts) <- self$featureData$FEATURE_ID

        distmat <- switch(
          metric,
          "unifrac" = rbiom::bdiv_distmat(
            biom = counts,
            bdiv = metric,
            weighted = weighted,
            tree = self$treeData,
            cpus = cpus
            ),
          "manhattan" = ,
          "euclidean" = ,
          "jaccard" = ,
          "bray" = rbiom::bdiv_distmat(
            biom = counts,
            bdiv = metric,
            weighted = weighted,
            cpus = cpus
            )
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
      # Switch case to compute relevant statistics
      stats_results <- switch(
        method,
        "pcoa" = pairwise_adonis(distmat, groups = self$metaData[[ group_by ]], perm = perm),
        "nmds" = pairwise_anosim(distmat, groups = self$metaData[[ group_by ]], perm = perm)
      )
      plot_list$anova_data <- stats_results

      # Normalization of eigenvalues
      if (method == "pcoa") {
        pcs$eig_norm <- pcs$eig %>%
          purrr::map(function(x) x / sum(pcs$eig) * 100) %>%
          unlist()

        # Collects loading scores into dataframe
        df_pcs_points <- data.table::data.table(pcs$points)
        colnames(df_pcs_points) <- paste0("PC", 1:ncol(df_pcs_points))
      } else if (method == "nmds") {
        df_pcs_points <- data.table::data.table(pcs$points)
        df_pcs_points$stress <- pcs$stress
      }
      plot_list$pcs <- df_pcs_points

      # Adds relevant data
      df_pcs_points[, groups := self$metaData[[ group_by ]] ]
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
          dist_explained = pcs$eig_norm,
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
    #' @param foldchange.threshold A numeric value used as a fold-change threshold to label and color significantly expressed features (Default: 0.06).
    #' @param abundance.threshold A numeric value used as an abundance threshold to size the scatter dots based on their mean relative abundance (default: 0.01).
    #' @param normalize A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @examples
    #' library(ggplot2)
    #' 
    #' obj_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' obj <- readRDS(obj_path)
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
                   foldchange.threshold = 0.06,
                   abundance.threshold = 0
                   ) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank) && length(feature_rank) != 1)
        cli::cli_abort("{feature_rank} needs to be a character with a length of 1")

      if (!is.character(condition.group) && length(condition.group) != 1) {
        cli::cli_abort("{condition.group} needs to be a character with a length of 1")
      } else if (!column_exists(condition.group, self$metaData)) {
        cli::cli_abort("{condition.group} does not exist in the metaData or is empty.")
      }
      if (!is.character(condition_A))
        cli::cli_abort("{condition_A} needs to be a character.")

      if (!is.character(condition_B))
        cli::cli_abort("{condition_B} needs to be a character.")

      if (!is.numeric(pvalue.threshold))
        cli::cli_abort("{pvalue.threshold} need to be numeric.")

      if (!is.numeric(foldchange.threshold))
        cli::cli_abort("{foldchange.threshold} need to be numeric.")

      if (paired && is.null(self$.samplepair_id)) {
        cli::cli_alert_warning("Paired is set to {paired} but SAMPLEPAIR_ID does not exist in the metaData.\n Differential feature analysis will not be done with paired set to FALSE!")
        paired <- FALSE
      }

      ## MAIN
      #--------------------------------------------------------------------#

      # Final output
      plot_list <- list()

      # Copies object to prevent modification of omics class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )
      
      # Restores omics class components
      on.exit(private$tmp_restore(), add = TRUE)

      # normalization if applicable
      if (normalize)
        self$normalize()

      # Subset by missing values
      self$removeNAs(condition.group)

      # Subset by samplepair completion
      if (paired && !is.null(self$.samplepair_id))
        self$samplepair_subset()

      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_merge(feature_rank = feature_rank,
                        feature_filter = feature_filter)

      # Extract mean relative abundance
      rel_abun <- as.matrix(Matrix::rowSums(self$countData) / ncol(self$countData))
      rownames(rel_abun) <- self$featureData[[ feature_rank ]]

      # Get data.table format relative abundances
      dt <- sparse_to_dtable(self$countData)[, (feature_rank) := self$featureData[[feature_rank]]]

      # Compute 2-fold expression based on (un)paired samples
      # Computes on equation of log2(A) - log2(B)
      # Supports multiple inputs for A and B.
      condition.labels <- data.table::setorderv(self$metaData,
                                                cols = c(self$.sample_id, condition.group))[[ condition.group ]]

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

        # Add relative abundance, and save data as output list
        dfe <- dfe[, "rel_abun" := rel_abun]
        plot_list$data <- dfe

        # Create & save volcano plot
        n_diff_columns <- sum(grepl("^Log2FC_", colnames(dfe)))

        plot_list$volcano_plot <- lapply(1:n_diff_columns, function(i) {
          volcano_plot(data = dfe,
                       logfold_col = paste0("Log2FC_", i),
                       pvalue_col = paste0("pvalue_", i),
                       feature_rank = feature_rank,
                       abundance_col = "rel_abun",
                       pvalue.threshold = pvalue.threshold,
                       logfold.threshold = foldchange.threshold,
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
    #' For now only works with headers that start with prefix `CONTRAST_`.
    #' @param feature_ranks A character vector as input to [`rankstat()`](#method-rankstat).
    #' @param feature_contrast A character vector of feature columns in the `featureData` to aggregate via [`feature_merge()`](#method-feature_merge).
    #' @param feature_filter A character vector to filter unwanted features, default: \code{c("uncultured")}
    #' @param distance_metrics A character vector specifying what (dis)similarity metrics to use, default \code{c("unifrac")}
    #' @param beta_div_table A path to pre-computed distance matrix, expects tsv/csv/txt file.
    #' @param alpha_div_table A path to pre-computed alpha diversity with rarefraction depth, expects tsv/csv/txt from qiime2, see \link{read_rarefraction_qiime}.
    #' @param normalize A boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @param weighted A boolean value, whether to compute weighted or unweighted dissimilarities (Default: TRUE).
    #' @param pvalue.threshold A numeric value, the p-value is used to include/exclude composition and foldchanges plots coming from alpha- and beta diversity analysis (Default: 0.05).
    #' @param perm A wholenumber, number of permutations to compare against the null hypothesis of \link[vegan]{adonis2} or \link[vegan]{anosim} (default: \code{perm=999}).
    #' @param cpus Number of cores to use, only used in [`ordination()`](#method-ordination) when beta_div_table is not supplied.
    #' @param filename A character to name the HTML report, it can also be a filepath (e.g. \code{"/path/to/report.html"}). Default: "report.html" in your current work directory.
    #' @importFrom patchwork plot_layout wrap_plots
    #' @return A report in HTML format
    autoFlow = function(feature_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                        feature_contrast = c("Phylum", "Family", "Genus"),
                        feature_filter = c("uncultured"),
                        distance_metrics = c("unifrac"),
                        beta_div_table = NULL,
                        alpha_div_table = NULL,
                        normalize = TRUE,
                        weighted = TRUE,
                        pvalue.threshold = 0.05,
                        perm = 999,
                        cpus = 1,
                        filename = paste0(getwd(), "/report.html")
                      ) {
    ## Error handling
    #--------------------------------------------------------------------#

    if (!is.character(filename) && length(filename) != 1)
      cli::cli_abort("{filename} needs to be a character with a length of 1")
      
    if (!is.character(feature_contrast) && length(feature_contrast) != 1) {
      cli::cli_abort("{feature_contrast} needs to be a character with a length of 1")
    } else if (!column_exists(feature_contrast, self$featureData)) {
      cli::cli_abort("{feature_contrast} does not exist in featureData!")
    }

    if (!is.null(beta_div_table) && !is.character(beta_div_table) && length(beta_div_table) != 1) {
      cli::cli_abort("{beta_div_table} needs to be a character with a length of 1")
    
      if (file.exists(beta_div_table))
        cli::cli_abort("{beta_div_table} already exists!")
    }

    if (!is.null(alpha_div_table) && !is.character(alpha_div_table) && length(alpha_div_table) != 1) {
      cli::cli_abort("{alpha_div_table} needs to be a character with a length of 1")

      if (file.exists(alpha_div_table))
        cli::cli_abort("{alpha_div_table} already exists!")
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
    private$tmp_link(
      .countData = self$countData,
      .featureData = self$featureData,
      .metaData = self$metaData,
      .treeData = self$treeData
    )
    
    # Restores omics class components
    on.exit(private$tmp_restore(), add = TRUE)

    # Collect columns: CONTRAST_ and VARIABLE_
    metacols <- colnames(self$metaData)

    CONTRAST_data <- self$metaData[, .SD, .SDcols = grepl("CONTRAST_", metacols)]
    CONTRAST_names <- colnames(CONTRAST_data)

    VARIABLE_data <- self$metaData[, .SD, .SDcols = grepl("VARIABLE_", metacols)]
    VARIABLE_names <- colnames(VARIABLE_data)

    #
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
    plots$rankstat_plot <- self$rankstat(feature_ranks)

    # Main loop
    if (CONTRAST_ncol > 0) {

      # Load custom distance matrix if supplied
      if (!is.null(beta_div_table)) {
        beta_div_table <- check_matrix(filepath = beta_div_table)
        beta_div_table <- beta_div_table[self$metaData[[self$.sample_id]], self$metaData[[self$.sample_id]]]
      }

      # Load custom rarefraction alpha diversity table if supplied
      if (!is.null(alpha_div_table)) {
        alpha_div_table <- read_rarefraction_qiime(filepath = alpha_div_table)
      }

      # Initialize plot containers
      composition_plots <- matrix(list(), CONTRAST_ncol, feature_nrow)
      Log2FC_plots <- matrix(list(), CONTRAST_ncol, feature_nrow)
      alpha_div_plots <- list()
      metrics_nrow <- length(distance_metrics)
      pcoa_plots <- matrix(list(), CONTRAST_ncol, metrics_nrow)
      nmds_plots <- matrix(list(), CONTRAST_ncol, metrics_nrow)

      # Initialize data containers
      composition_data <- matrix(list(), CONTRAST_ncol, feature_nrow)
      Log2FC_data <- matrix(list(), CONTRAST_ncol, feature_nrow)
      alpha_div_data <- list()
      pcoa_data <- matrix(list(), CONTRAST_ncol, metrics_nrow)
      nmds_data <- matrix(list(), CONTRAST_ncol, metrics_nrow)

      for (i in 1:CONTRAST_ncol) {
        col_name <- CONTRAST_names[i]
        conditions <- NULL
        cli::cli_alert_info(paste0("Processing ... column: ", col_name, " \n"))

        #--------------------------------------------------------------------#
        ## Alpha diversity
        #--------------------------------------------------------------------#
        if (inherits(alpha_div_table, "data.table")) {
          dt_final <- base::merge(alpha_div_table,
                                  self$metaData[, .SD, .SDcols = c(self$.sample_id, col_name)],
                                  by = self$.sample_id,
                                  all.x = TRUE) %>%
            na.omit(cols = col_name)

          res <- diversity_plot(
            data = dt_final,
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
                paired = ifelse(!is.null(self$.samplepair_id), TRUE, FALSE)
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
              cpus = cpus
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
              perm = perm
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
          nmds_plots[[i, j]] <- patchwork::wrap_plots(res[c("anova_plot", "scores_plot")],
                                                      nrow = 1) +
            patchwork::plot_layout(widths = c(rep(5, 3)),
                                   guides = "collect")
          nmds_data[[i, j]] <- list(
            stats = res$anova_data,
            dist_mat = res$dist,
            pcs = res$pcs
          )
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
                paired = ifelse(!is.null(self$.samplepair_id), TRUE, FALSE),
                normalize = normalize,
                condition.group = col_name,
                condition_A = c(conditions$group1),
                condition_B = c(conditions$group2)
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
                  condition_B = c(conditions$group2)
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
      plots$nmds_plots <- is_empty(nmds_plots)

      # Checks if data aren't empty
      data$composition_data <- is_empty(composition_data)
      data$Log2FC_data <- is_empty(Log2FC_data)
      data$alpha_div_data <- is_empty(alpha_div_data)
      data$pcoa_data <- is_empty(pcoa_data)
      data$nmds_data <- is_empty(nmds_data)
    }
    
    #--------------------------------------------------------------------#
    ## CREATING REPORT
    #--------------------------------------------------------------------#

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
  }
  ),
  private = list(
    # Creates a temporary save of self components
    tmp_store = NULL,
    tmp_link = function(.countData = NULL, .featureData = NULL, .metaData = NULL, .treeData = NULL) {
      private$tmp_store <<- list(
                            .countData = .countData,
                            .metaData = .metaData,
                            .featureData = .featureData,
                            .treeData = .treeData
                            )
    },
    tmp_restore = function() {
      # Restores self components if applicable!
      if (!is.null(private$tmp_store$.countData)) self$countData <- private$tmp_store$.countData
      if (!is.null(private$tmp_store$.metaData)) self$metaData <- private$tmp_store$.metaData
      if (!is.null(private$tmp_store$.featureData)) self$featureData <- private$tmp_store$.featureData
      if (!is.null(private$tmp_store$.treeData)) self$treeData <- private$tmp_store$.treeData
      return(invisible(self))
    },
    check_table = function(data) {
    #
    # Creates a data.table from given filepath, data.frame or data.table
    #
    if (is.character(data) && length(data) == 1 && file.exists(data))
      return(data.table::fread(data, header = TRUE))

    if (inherits(data, "data.table"))
      return(data)

    if (is.data.frame(data))
      return(data.table::as.data.table(data))

    stop("Input must be a filepath, data.frame, or data.table.")
  },
  check_matrix = function(data) {
    #
    # Creates a sparseMatrix from given filepath, matrix or sparseMatrix
    #
    if (is.character(data) && length(data) == 1 && file.exists(data)) {
      # 
      # If filepath originates from an excel file, it may contain trailing spaces, or letters, which are removed.
      #
      dt <- data.table::fread(data, header = TRUE)
      # Change character values to numeric
      for (col in names(dt)) {
        dt[is.na(get(col)), (col) := 0]
        dt[get(col) == "", (col) := 0]
      }

      # Convert to matrix format
      mat <- Matrix::Matrix(
        data = as.matrix(dt),
        dimnames = list(rownames(dt), colnames(dt))
      )
      
      # Return sparseMatrix
      return(as(mat, "sparseMatrix"))
    }

    if (inherits(data, "sparseMatrix"))
      return(data)

    if (is.matrix(data))
      sp_mat <- Matrix::Matrix(data, sparse = TRUE)
      return(sp_mat)

    stop("Input must be a filepath, matrix, or sparseMatrix.")
    }
  )
)
