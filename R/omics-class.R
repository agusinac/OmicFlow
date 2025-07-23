#' Abstract 'omics' class
#'
#' @description This is the abstract class 'omics', contains a variety of methods that are inherited and applied in the omics classes:
#' \link[OmicFlow]{metagenomics}, proteomics and metabolomics.
#'
#' @details
#' Every class is created with the \link[R6]{R6Class} method. Methods are either public or private, and only the public components are inherited by other omics classes.
#' The omics class by default uses triplet \link[data.table]{data.table} data structures for quick and efficient data manipulation and returns the object by reference, same as the R6 class.
#' The method by reference is very efficient when dealing with big data.
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

    #' @field .feature_id String value, default name for the feature identifiers.
    .feature_id = "FEATURE_ID",

    #' @field .sample_id String value, default name for the sample identifiers.
    .sample_id = "SAMPLE_ID",

    #' @field .samplepair_id String value, default name for the sample pair identifiers.
    .samplepair_id = "SAMPLEPAIR_ID",

    #' @description
    #' Wrapper function that is inherited and adapted for each omics class.
    #' The omics classes requires a metadata samplesheet, that is validated by the metadata_schema.json. It requires a column `SAMPLE_ID` and optionally a `SAMPLEPAIR_ID` or `FEATURE_ID` can be supplied. The `SAMPLE_ID` will be used to link the metaData to the countData, and will act as the key during subsetting of other columns.
    #' To create a new object use \code{omics$new()}
    #' @param countData countData A path to an existing file, data.table or data.frame.
    #' @param featureData A path to an existing file, data.table or data.frame.
    #' @param metaData A path to an existing file, data.table or data.frame.
    #' @return A new `omics` object.
    #'
    initialize = function(countData = NA, featureData = NA, metaData = NA) {
      #-------------------#
      ###   metaData    ###
      #-------------------#
      if (!is.na(metaData)) {
        self$metaData <- data.table::fread(metaData, header = TRUE)
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
          if (column_exists(self$.samplepair_id, self$metaData))
            self$.samplepair_id <- NULL

        } else {
          errors <- attr(self$.valid_schema, "errors")
          cli::cli_abort(
            "JSON validation failed: \n{ paste(errors$message, collapse = '\n')}"
            )
        }

      } else {
        cli::cli_abort(c(
          "metaData cannot be empty, please provide a tab or comma separated file"
        ))
      }

      #-------------------#
      ###  featureData  ###
      #-------------------#
      if (!is.na(featureData)) {
        self$featureData <- data.table::fread(featureData,
                                              header = TRUE)

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
      if (!is.na(countData)) {
        self$countData <- read_sparseTable(countData)
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
    #' Validates an input metadata against the JSON schema, see [`metadata_schema.json`](./metadata_schema.json).
    #' This function is used when a metaData table is supplied and automatically checked during the creation of a new object.
    #'
    validate = function() {
      # Creates temporary json file from metadata
      tmp_json <- base::tempfile(fileext = ".json")

      json_data <- jsonlite::toJSON(self$metaData,
                                    dataframe = "rows",
                                    pretty = TRUE,
                                    auto_unbox = TRUE)

      writeLines(json_data, tmp_json)

      self$.valid_schema <- jsonvalidate::json_validate(
        tmp_json, "metadata_schema.json",
        engine = "ajv",
        verbose = TRUE,
        error = FALSE,
        strict = TRUE
        )

      unlink(tmp_json)
      invisible(self)
    },
    #' @description
    #' Removes empty (zero) values by row and column from the `countData`
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$removeZeros()
    #'
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
    #' Remove NAs from metaData and updates the countData object fields.
    #' @param column The column from where NAs should be removed, this can be either integers or strings.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$removeNAs(column = "treatment")
    #'
    removeNAs = function(column) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (is.integer(column) && length(column) <= length(colnames(self$metaData))) {
        column <- colnames(self$metaData[column])

      } else {
        cli::cli_abort("{column} indexes are out of bound")
      }

      if (!column_exists(column, self$metaData))
        cli::cli_abort("{column} do not exist in the metaData or one of the specified columns is completely empty!")

      ## MAIN
      #--------------------------------------------------------------------#

      self$metaData <- na.omit(self$metaData, cols = column)
      self$countData <- self$countData[, self$metaData[[ self$.sample_id ]]]
      invisible(self)
    },
    #' @description
    #' Feature subset (based on featureData), automatically applies \code{removeZeros}
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in `featureData`.
    #' Only rows for which all conditions evaluate to TRUE are kept.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$feature_subset(rank1 == "Streptococcus")
    #' obj$feature_subset(rank1 %in% c("Streptococcus", "uncultured"))
    #'
    feature_subset = function(...) {
      rows_to_keep <- self$featureData[, ...]
      self$featureData <- self$featureData[rows_to_keep, ]
      self$countData <- self$countData[rows_to_keep, ]
      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Sample subset (based on metaData), automatically applies \code{removeZeros}
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in `metaData`.
    #' Only rows for which all conditions evaluate to TRUE are kept.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$sample_subset(cycle == "t1")
    #' obj$sample_subset(cycle %in% c("t1", "t5"))
    #'
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
    #' Samplepair subset (based on metaData), automatically applies \code{removeZeros}
    #' @param num_unique_pairs An integer value to define the number of pairs to subset. The default is NULL, meaning the maximum number of unique pairs will be used to subset the data. Let's say you have three samples for each pair, then the `num_unique_pairs` will be set to 3.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$samplepair_subset()
    #' obj$samplepair_subset(num_unique_pairs = 2)
    #'
    samplepair_subset = function(num_unique_pairs = NULL) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.null(num_unique_pairs) && !is.integer(num_unique_pairs))
        cli::cli_abort("{num_unique_pairs} must contain integers!")

      ## MAIN
      #--------------------------------------------------------------------#

      if (is.null(num_unique_pairs)) {
        counts <- self$metaData[, .(unique_count = data.table::uniqueN(SAMPLE_ID)), by = SAMPLEPAIR_ID]
        num_unique_pairs <- counts[, max(unique_count)]
      }

      self$metaData <- self$metaData[SAMPLEPAIR_ID %in% counts[unique_count == num_unique_pairs, SAMPLEPAIR_ID]]
      self$countData <- self$countData[, self$metaData[[self$.sample_id]] ]
      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Agglomerates features by column, automatically applies \code{removeZeros}.
    #' @param feature_rank A character value or vector of columns to aggregate from the `featureData`.
    #' @param feature_filter A character value or vector of characters to remove features via regex pattern.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$feature_glom(feature_rank = c("Kingdom", "Phylum"))
    #' obj$feature_glom(feature_rank = "Genus", feature_filter = c("uncultured", "metagenome"))
    #'
    feature_glom = function(feature_rank, feature_filter = NA) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank))
        cli::cli_abort("{feature_rank} needs to be a character or vector containing characters")

      if (!is.na(feature_filter) && !is.character(feature_filter))
        cli::cli_abort("{feature_filter} needs to be a character or vector containing characters")

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


      if (!is.na(feature_filter)) {
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
    #' Performs transformation on countData as a Triplet sparse matrix \link[Matrix]{uniqTsparse}
    #' @param FUN A function such as \code{log2}, \code{log}
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$transform(log2)
    #'
    transform = function(FUN) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!inherits(FUN, "function"))
        cli::cli_abort("{FUN} must be a function!")

      ## MAIN
      #--------------------------------------------------------------------#

      self$countData@x <- fun(self$countData@x)
      invisible(self)
    },
    #' @description
    #' Relative abundance computation by column sums.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$normalize()
    #'
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
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' plt <- obj$rankstat()
    #' plt
    #' @return A \link[ggplot2]{ggplot} object.
    #'
    rankstat = function(feature_ranks) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (is.integer(feature_ranks) && length(feature_ranks) <= length(colnames(self$featureData))) {
        column <- colnames(self$featureData[feature_ranks])

      } else {
        cli::cli_abort("{column} indexes are out of bound.")
      }

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
                    y = "Number of ASVs classified"))
    },
    #' @description
    #' Alpha diversity based on \link[OmicFlow]{diversity}
    #' @param col_name The metaData column of categorical variables to create a ggplot object.
    #' @param method Diversity metric such as "shannon", "invsimpson" or "simpson"
    #' @param Brewer.palID Palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link[OmicFlow]{fetch_palette}.
    #' @param evenness A boolean wether to divide diversity by number of species, see \link[vegan]{specnumber}.
    #' @param paired A boolean value to perform paired analysis in \link[stats]{wilcox.test} and samplepair subsetting via \link[OmicFlow]{samplepair_subset}
    #' @param p.adjust.method A character variable to specify the p.adjust.method to be used, default is 'fdr'.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' plt <- obj$alpha_diversity(col_name = "treatment",
    #'                            method = "shannon")
    #'
    #' @return A \link[ggplot2]{ggplot} object.
    #' @seealso \link[OmicFlow]{diversity_plot}
    alpha_diversity = function(col_name,
                               method = c("shannon", "invsimpson", "simpson"),
                               Brewer.palID = "Set2",
                               evenness = FALSE,
                               paired = FALSE,
                               p.adjust.method = "fdr") {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(col_name) && length(col_name) != 1) {
        cli::cli_abort("{col_name} must be a character and of length 1")
      } else if (!column_exists(column_name, self$metaData)) {
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

      # Subset by samplepair completion
      if ( paired && !is.null(self$.samplepair_id) )
        self$samplepair_subset()

      # Alpha diversity based on 'method'
      div <- data.table::data.table(diversity(x = self$countData, index=method))
      div[, (col_name) := self$metaData[, .SD, .SDcols = c(col_name)]]
      # Adjusts for evenness
      if (evenness) div$V1 <- div$V1 / log(vegan::specnumber(div$V1))

      # get colors
      colors <- fetch_palette(self$metaData, col_name, Brewer.palID)

      # Create and saves plots
      plot_list$data <- div
      plot_list$diversity <- diversity_plot(data = na.omit(div),
                                            values = "V1",
                                            col_name = col_name,
                                            palette = colors,
                                            method = method,
                                            paired = paired,
                                            p.adjust.method = p.adjust.method)

      # Restores omics class components
      private$tmp_restore()

      return(plot_list)
    },
    #' @description
    #' Visualization of compositional data.
    #' @param feature_rank A featureData column name to visualize.
    #' @param feature_filter Removes features by name, works on single strings or vector of strings.
    #' @param col_name Optional, a string or vector of strings to add to the final compositional data output.
    #' @param feature_top Integer of the top features to visualize, the max is 15, due to a limit of palettes.
    #' @param Brewer.palID Palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link[OmicFlow]{fetch_palette}.
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #'
    #' result <- obj$composition(feature_rank = "Genus",
    #'                           feature_filter = c("uncultured"),
    #'                           feature_top = 10)
    #'
    #' plt <- composition_plot(data = result$data,
    #'                         palette = result$palette,
    #'                         feature_rank = "Genus")
    #'
    #' @return A long \link[data.table]{data.table} table.
    #' @seealso \link[OmicFlow]{composition_plot}
    composition = function(feature_rank,
                           feature_filter = NA,
                           col_name = NULL,
                           feature_top = c(10, 15),
                           Brewer.palID = "RdYlBu",
                           remove_na = FALSE) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.null(col_name) && !is.character(col_name) && length(col_name) != 1) {
        cli::cli_abort("{col_name} must be a character and of length 1")
      } else if (!column_exists(col_name, self$metaData)) {
        cli::cli_abort("The specified {col_name} does not exist in the metaData.")
      }

      if (!is.integer(feature_top)) {
        cli::cli_abort("{feature_top} must be an integer!")
      } else if (feature_top > 15) {
        cli::cli_alert_warning("The {feature_top} is set to an integer higher than 15.\n This may lead that colors are difficult to be distinguished.\n For color-blind people it is recommended to use a feature_top of maximum 15.")
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

      # Agglomerate by feature_rank
      self$feature_glom(feature_rank = feature_rank, feature_filter = feature_filter)

      # Normalizes sample counts
      self$normalize()

      # Remove NAs when col_name is specified
      if (!is.null(col_name) & remove_na)
        self$removeNAs(col_name)

      # Convert sparse matrix to data.table (safe since feature_glom shrinks the sparse matrix)
      counts <- data.table::data.table(as.matrix(self$countData))

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
        chosen_palette <- viridis::viridis(df_taxa_len - 1)
      } else if (df_taxa_len-1 <= 15 & df_taxa_len-1 > 10) {
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

      # Restores omics class components
      private$tmp_restore()

      # returns results as list
      return(
        list(
          data = composition_final,
          palette = taxa_colors_ordered
        )
      )
    },
    #' @description
    #' Ordination of countData with statistical tests.
    #' @param metric A dissimilarity or similarity metric to be applied on the `countData`, thus far supports 'bray', 'jaccard' and 'unifrac' when a tree is provided via `treeData`.
    #' @param method Ordination method, supports "pcoa" and "nmds".
    #' @param distmat A custom distance matrix in either \link[stats]{dist} or \link[Matrix]{Matrix} format.
    #' @param group_by A string variable as metaData column to be used for the PERMANOVA or ANOSIM statistical test.
    #' @param weighted A Boolean value, whether to compute weighted or unweighted dissimilarities (Default: TRUE).
    #' @param normalize A Boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @param cpus An Integer value, indicating the number of processes to spawn (Default: 1).
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #'
    #' pcoa_plots <- obj$ordination(metric = "bray",
    #'                              method = "pcoa",
    #'                              group_by = "treatment",
    #'                              weighted = TRUE,
    #'                              parallel = TRUE,
    #'                              normalize = TRUE)
    #' pcoa_plots
    #'
    #' @return A list of \link[ggplot2]{ggplot} object.
    #' @seealso \link[OmicFlow]{ordination_plot}, \link[OmicFlow]{stats_plot}, \link[OmicFlow]{pairwise_anosim}, \link[OmicFlow]{pairwise_adonis}
    ordination = function(metric = c("bray", "jaccard", "unifrac"),
                          method = c("pcoa", "nmds"),
                          group_by,
                          distmat = NULL,
                          weighted = TRUE,
                          normalize = TRUE,
                          cpus = 1,
                          perm=999) {

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

      if (!is.integer(cpus))
        cli::cli_abort("{cpus} need to be an integer!")

      if (!is.integer(perm))
        cli::cli_abort("Permutations {perm} need to be an integer")

      if (!inherits(distmat, "Matrix") || !inherits(distmat, "dist"))
        cli::cli_abort("custom distance matrix (distmat) need to be of class Matrix or dist")

      ## MAIN
      #--------------------------------------------------------------------#

      # Copies object to prevent modification of omics class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

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
      counts <- slam::as.simple_triplet_matrix(self$countData)
      rownames(counts) <- self$featureData$FEATURE_ID

      distmat <- switch(
        metric,
        "unifrac" = rbiom::bdiv_distmat(biom = counts,
                                        bdiv = metric,
                                        weighted = weighted,
                                        tree = self$treeData,
                                        cpus = cpus),
        "manhattan" = ,
        "euclidean" = ,
        "jaccard" = ,
        "bray" = rbiom::bdiv_distmat(biom = counts,
                                     bdiv = metric,
                                     weighted = weighted,
                                     cpus = cpus)
      )

      plot_list$distances <- distmat


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
      plot_list$pcs <- pcs

      # Adds relevant data
      df_pcs_points[, groups := self$metaData[[ group_by ]] ]
      df_pcs_points[, samples := row.names(df_pcs_points) ]

      # # Pairwise dimensions
      # if (pca.pairwise & method == "pcoa") {
      #   # Finds number of dimensions that explain 80% of distances
      #   n_dimensions = 0
      #   sum_eig = 0
      #   for (eig in pcs$eig_norm) {
      #     if (sum_eig < pca.max.explained) {
      #       sum_eig <- sum_eig + eig
      #       n_dimensions <- n_dimensions + 1
      #     } else break
      #   }
      #
      #   # Creates paired combinations of dimensions into a list of plots
      #   n_dim_pairs <- utils::combn(seq(n_dimensions), 2)
      #   pdf(paste0(outdir, "/pairwise_PCoA.pdf"))
      #   for (i in seq(ncol(n_dim_pairs))) {
      #     pair <- n_dim_pairs[, i]
      #     print(pcoa_plot(df_pcs_points, pcs, pair, metric))
      #   }
      #   dev.off()
      # }

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
          labs(title = "Screeplot of first 10 PCs",
               x = "Principal Components (PCs)",
               y = "dissimilarity explained [%]")

        # PERMANOVA
        plot_list$anova_plot <- stats_plot(stats_results,
                                           X = "pairs",
                                           Y = "F.Model",
                                           Label = "p.adj",
                                           Y_title = "Pseudo F test statistic",
                                           plot.title = "PERMANOVA")
        # Loading score plot
        plot_list$scores_plot <- ordination_plot(df_pcs_points,
                                                 pcs,
                                                 pair=c("PC1", "PC2"),
                                                 metric)

      } else if (method == "nmds") {
        plot_list$anova_plot <- stats_plot(stats_results,
                                           X = "pairs",
                                           Y = "anosimR",
                                           Label = "p.adj",
                                           Y_title = "ANOSIM R statistic",
                                           plot.title = "ANOSIM")

        plot_list$scores_plot <- ordination_plot(df_pcs_points,
                                                 pcs,
                                                 pair=c("MDS1", "MDS2"),
                                                 metric)
      }

      # Restores omics class components
      private$tmp_restore()

      return(plot_list)
    },
    #' @description
    #' Differential feature expression using the \link[OmicFlow]{foldchange} for both paired and non-paired samples.
    #' @param feature_rank A character value or vector of columns to aggregate from the `featureData`.
    #' @param feature_filter A character value or vector of characters to remove features via regex pattern (Default: NULL).
    #' @param feature_top Integer of the top features to visualize (Default: NULL, everything will be used).
    #' @param paired A Boolean value, the paired is only applicable when a `SAMPLEPAIR_ID` column exists within the `metaData`. See \link[stats]{wilcox.test}
    #' @param condition.group A string value for an existing column name in `metaData`, wherein the conditions A and B are located.
    #' @param condition_A A character value or vector of characters.
    #' @param condition_B A character value or vector of characters.
    #' @param pvalue.threshold An Integer value, a P-value threshold to label and color significant features (Default: 0.05).
    #' @param foldchange.threshold An Integer value, a fold-change threshold to label and color significantly expressed features (Default: 0.06).
    #' @param normalize A Boolean value, whether to [`normalize()`](#method-normalize) by total sample sums (Default: TRUE).
    #' @examples
    #' obj <- omics$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #'
    #' unpaired <- obj$differential_feature_expression(feature_rank = "Genus",
    #'                                                 paired = FALSE,
    #'                                                 condition.group = "treatment",
    #'                                                 condition_A = c("H"),
    #'                                                 condition_B = c("T"))
    #'
    #' paired <- obj$differential_feature_expression(feature_rank = "Genus",
    #'                                               paired = TRUE,
    #'                                               condition.group = "cycle",
    #'                                               condition_A = c("t2", "t3"),
    #'                                               condition_B = c("t1", "t2"),
    #'                                               feature_top = 20)
    #'
    #' @return
    #' * A list of \link[ggplot2]{ggplot} object.
    #' * A long \link[data.table]{data.table} table.
    #' @seealso \link[OmicFlow]{volcano_plot}, \link[OmicFlow]{ViolinBoxPlot}, \link[OmicFlow]{paired_fold}, \link[OmicFlow]{unpaired_fold}
    differential_feature_expression = function(feature_rank,
                                               feature_filter = NULL,
                                               feature_top = NULL,
                                               paired = FALSE,
                                               normalize = TRUE,
                                               condition.group,
                                               condition_A,
                                               condition_B,
                                               pvalue.threshold = 0.05,
                                               foldchange.threshold = 0.06
                                               ) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank) && length(feature_rank) != 1)
        cli::cli_abort("{feature_rank} needs to be a character with a lenght of 1")

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

      # Subset by samplepair completion
      if ( paired && !is.null(self$.samplepair_id) )
        self$samplepair_subset()

      # Subset by missing values
      self$removeNAs(condition.group)

      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_glom(feature_rank = feature_rank,
                        feature_filter = feature_filter)

      # normalization if applicable
      if (normalize)
        self$normalize()

      # Check how many features to select (depended if volcano is desired)
      if (!is.null(feature_top)) {
        feature_top <- feature_top
      } else {
        feature_top <- nrow(self$featureData)
      }

      # Extract relative abundance
      rel_abun <- as.matrix(Matrix::rowMeans(self$countData[1:feature_top,]))
      rownames(rel_abun) <- self$featureData[[ feature_rank ]]

      # Creates long table of relative abundance
      dt <- sparse_to_dtable(self$countData)[, (feature_rank) := self$featureData[[feature_rank]]]
      stats_dt <- base::merge(data.table::melt(dt,
                                               measure.vars = colnames(dt)[!grepl(feature_rank, colnames(dt))],
                                               variable.name = self$.sample_id,
                                               value.name = "values"),
                              self$metaData[, .SD, .SDcols = c(self$.sample_id, condition.group)],
                              by = self$.sample_id)

      # Create row_sums
      dt[, row_sum := rowSums(.SD), .SDcols = !c(feature_rank)]

      # Orders by row_sum in descending order
      countTable <- data.table::setorder(dt, -row_sum)[1:feature_top, .SD, .SDcols = !c("row_sum")]
      features <- countTable[[ feature_rank ]]
      self$countData <- as(as.matrix(countTable[, .SD, .SDcols = !c(feature_rank)]), "sparseMatrix")

      # Log2 transform taxa
      self$transform(log2)

      # Subset by top features
      stats_dt <- stats_dt[stats_dt[[feature_rank]] %in% features]
      dt <- sparse_to_dtable(self$countData)[, (feature_rank) := features]

      # Compute 2-fold expression based on (un)paired samples
      # Computes on equation oflog2(A) - log2(B)
      # Supports multiple inputs for A and B.
      condition.labels <- data.table::setorderv(self$metaData,
                                                cols = c(self$.sample_id, condition.group))[[ condition.group ]]

      # paired samples
      DFE <- foldchange(
        data = dt,
        sample.id = self$.sample_id,
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
        DFE <- DFE[, "rel_abun" := rel_abun]
        plot_list$data <- DFE

        # Create & save volcano plot
        n_diff_columns <- sum(grepl("^Log2FC_", colnames(DFE)))

        plot_list$volcano_plot <- lapply(1:n_diff_columns, function(i) {
          volcano_plot(data = DFE,
                       X = paste0("Log2FC_", i),
                       Y = paste0("pvalue_", i),
                       feature_rank = feature_rank,
                       pvalue.threshold = pvalue.threshold,
                       logfold.threshold = foldchange.threshold,
                       label_A = condition_A,
                       label_B = condition_B)
        })

      # Restores omics class components
      private$tmp_restore()

      return(plot_list)
    },
    # triplot = function(feature_rank,
    #                    feature_filter = NA,
    #                    sample.id = self$.sample_id,
    #                    metadata.col = NA,
    #                    choice_dim = c("RDA1", "PC1"),
    #                    pairwise = FALSE,
    #                    Brewer.palID = "Set2",
    #                    counts.scalar = 1) {
    #   # Copies object to prevent modification of omics class components
    #   private$tmp_link(
    #     .countData = self$countData,
    #     .featureData = self$featureData,
    #     .metaData = self$metaData,
    #     .treeData = self$treeData
    #   )
    #
    #   ## Return list of components
    #   results <- list()
    #
    #   # Nested functions, later to be combined within omics-class
    #   logn <- function(otu_tab, scalar=1) {
    #     # log-transform, center
    #     # Y' = log ( A * Y + 1 ) ; where A is the 'strength' of the log transformation : 1, 10, 100, 1000, etc., default = 1
    #     otu_tab.log <- ( scalar * otu_tab ) + 1
    #     otu_tab.log <- log( otu_tab.log )
    #     otu_tab.sc <- scale(otu_tab.log, center = TRUE, scale = FALSE)
    #     return(otu_tab.sc)
    #   }
    #
    #   eigen_80 <- function(eig_explained) {
    #     sum_variance = 0
    #     counter = 1
    #     for (i in 1:length(eig_explained)) {
    #       sum_variance <- sum_variance + eig_explained[i]
    #       counter <- counter + 1
    #       if (sum_variance >= 80) break
    #     }
    #
    #     return(counter)
    #   }
    #
    #   subset_by_dimensions <- function(model, dimensions) {
    #     perc_explained <- round(100*(summary(model)$cont$importance[2, dimensions]),2)
    #     n_dim_pairs <- dimensions[1:eigen_80(perc_explained)]
    #     return(perc_explained)
    #   }
    #
    #   subset_by_species <- function(model, scores_species, pc) {
    #     species_explained <- utils::head(base::sort(round(100*scores_species[, pc]^2, 3), decreasing = TRUE))
    #     scores_species_explained <- scores_species[rownames(scores_species) %in% names(species_explained),]
    #
    #     result <- list(
    #       scores = scores_species_explained,
    #       explained_PC1 = species_explained
    #     )
    #
    #     return(result)
    #   }
    #
    #   # Main pairwise code
    #   if (pairwise == TRUE) {
    #     # Creates a vector of 10 dimensions (PC1 - PC10)
    #     pairwise_dims <- sprintf("PC%d", seq(1:11))
    #     subset.result <- subset_by_dimensions(model, pairwise_dims)
    #
    #     # save pdf
    #     pdf(paste0(outdir, "/pairwise_PCA.pdf"))
    #     for (i in seq(ncol(subset.result$n_dim_pairs))) {
    #       pair <- subset.result$n_dim_pairs[, i]
    #       scores_species_explained <- subset_by_species(model, scores_species, pc = pair[1])
    #
    #       triplot(model, target_col, self$metaData, subset.result$var_explained, scores_species,
    #               scores_species_explained, scores_sites,
    #               pc1 = pair[1],
    #               pc2 = pair[2])
    #     }
    #     dev.off()
    #   }
    #
    #   # Agglomerate taxa by feature rank and filter unwanted taxa
    #   self$feature_glom(feature_rank = feature_rank, feature_filter = feature_filter)
    #   self$normalize()
    #
    #   # Remove NAs
    #   if (!is.na(metadata.col)) {
    #     self$removeNAs(metadata.col)
    #   } else stop("metadata.col is empty!")
    #
    #   counts <- t(as.matrix(self$countData[, self$metaData[[ sample.id ]] ]))
    #   dimnames(counts)[[2]] <- self$featureData[[ feature_rank ]]
    #
    #   # Subsets user specified dimensions
    #   pc1 <- choice_dim[1]
    #   pc2 <- choice_dim[2]
    #
    #   # Transformation of counts and modelling to RDA
    #   counts.log <- logn(counts, scalar = counts.scalar)
    #   model <- vegan::rda(counts.log ~ get(metadata.col, self$metaData) + Condition(NULL),
    #                       data = self$metaData,
    #                       scale = FALSE,
    #                       na.action = na.fail,
    #                       subset = NULL)
    #
    #   results$model <- model
    #
    #   # MAIN code
    #   # Subset species and sites scores
    #   model.dims <- dim(model$CCA$u)[2] + dim(model$CA$u)[2]
    #   scores_species <- vegan::scores(x = model, display = "species", choices = c(1:model.dims), scaling=0)
    #   scores_sites <- vegan::scores(x = model, display = "sites", choices = c(1:model.dims))
    #
    #   # Subset species most fitted/captured by user defined dimensions
    #   choice_dim.scores_species_explained <- subset_by_species(model, scores_species, pc = choice_dim[1])
    #   choice_dim.explained <- subset_by_dimensions(model, choice_dim)
    #
    #   # Include relative abundance and significant groups in scores_sites
    #   rel_abun <- colSums(counts)
    #   Explained_species <- rownames(choice_dim.scores_species_explained$scores)
    #   scores_species_merged <- data.table::data.table(cbind(scores_species, rel_abun))
    #   scores_species_merged$taxa <- rownames(scores_species)
    #
    #   # Creating color palette
    #   chosen_palette <- RColorBrewer::brewer.pal(length(Explained_species), Brewer.palID)
    #   colors <- stats::setNames(chosen_palette, Explained_species)
    #
    #   # include groups for labelling and size
    #   scores_species_merged[, explained_species_label := ifelse(taxa %in% Explained_species, taxa, "")]
    #   scores_species_merged[, explained_species_size := ifelse(taxa %in% Explained_species, rel_abun, 0)]
    #
    #   #Fetch groups
    #   mygroups <- get(metadata.col, self$metaData)
    #   fills <- stats::setNames(RColorBrewer::brewer.pal(length(unique(mygroups)), "Set1"), unique(mygroups))
    #
    #   # to be named: scores_sites
    #   dt <- data.table::data.table(data.frame(pc1 = scores_sites[, pc1],
    #                                           pc2 = scores_sites[, pc2],
    #                                           group = mygroups))
    #
    #   results$data <- dt
    #   # Get centroid centers for annotation
    #   df_mean.ord <- stats::aggregate(dt, by=list(dt$group),mean)
    #   colnames(df_mean.ord) <- c("Group", "x", "y")
    #   df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% mygroups, ]
    #
    #   rslt.hull <- vegan::ordihull(scores_sites[, c(pc1, pc2)],
    #                                groups = mygroups,
    #                                show.group = mygroups,
    #                                draw = 'none')
    #
    #   # Initialize an empty data.table
    #   df_hull <- data.table::data.table(Group = character(), x = numeric(), y = numeric())
    #
    #   # Loop through groups and bind data
    #   for (g in mygroups) {
    #     g <- as.character(g)
    #     x <- rslt.hull[[g]][ , 1]
    #     y <- rslt.hull[[g]][ , 2]
    #     Group <- rep(g, length(x))
    #     df_temp <- data.table::data.table(Group = Group, x = x, y = y)
    #     df_hull <- rbind(df_hull, df_temp, use.names = TRUE, fill = TRUE)
    #   }
    #
    #   # Convert to data.table
    #   data.table::setDT(df_hull)
    #
    #   # Restores omics class components
    #   private$tmp_restore()
    #
    #
    #   results$plot <- ggplot() +
    #       # Polygon layer with first fill scale
    #       geom_polygon(data = df_hull,
    #                    aes(x = x,
    #                        y = y,
    #                        fill = Group),
    #                    alpha = 0.2,
    #                    color = "gray40",
    #                    show.legend = FALSE) +
    #       scale_fill_manual(values = fills) +
    #       ggrepel::geom_label_repel(data=df_mean.ord,
    #                                 aes(x=x, y=y, label=Group, fill=Group),
    #                                 color = "black",
    #                                 show.legend = FALSE) +
    #       guides(fill = "none") +
    #
    #       # Reset fill scale for points
    #       ggnewscale::new_scale_fill() +
    #
    #       # Main points layer with second fill scale
    #       geom_point(data = dt,
    #                  aes(x = .data[["pc1"]],
    #                      y = .data[["pc2"]]),
    #                  fill = "steelblue",
    #                  col = "black",
    #                  shape = 21) +
    #
    #       # Species points layer
    #       geom_point(data = scores_species_merged,
    #                  aes(x = .data[[pc1]],
    #                      y = .data[[pc2]],
    #                      size = .data[["explained_species_size"]],
    #                      col = .data[["explained_species_label"]],
    #                      stroke = ifelse(scores_species_merged$explained_species_label != "", 1.5, 0.5)),
    #                  show.legend = TRUE,
    #                  shape = 22) +
    #
    #       # Remaining formatting
    #       theme_bw() +
    #       theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    #             text = element_text(size = 12),
    #             legend.text = element_text(size = 12),
    #             legend.title = element_text(size = 14)) +
    #       scale_size_continuous(name = "rel. Abun.") +
    #       scale_color_manual(name = paste0(pc1, " explained species"),
    #                          values = colors) +
    #       labs(x = paste0(pc1, " (", choice_dim.explained[pc1], "%)"),
    #            y = paste0(pc2, " (", choice_dim.explained[pc2], "%)")) +
    #       guides(fill = guide_legend(position = "bottom", override.aes = list(size = 2, color = "white")),
    #              colour = guide_legend(override.aes = list(stroke = 1.5)))
    #
    #     return(results)
    # },
    # correlation = function(feature_rank, feature_filter = NA, sample.id = "SAMPLE-ID",
    #                        cor_method = "spearman", cor_columns = c("BMI", "Weight"), cor_threshold = 0.6, normalize = TRUE) {
    #   # Copies object to prevent modification of omics class components
    #   private$tmp_link(
    #     .countData = self$countData,
    #     .featureData = self$featureData,
    #     .metaData = self$metaData,
    #     .treeData = self$treeData
    #   )
    #
    #   # Agglomerate taxa by feature rank and filter unwanted taxa
    #   self$feature_glom(feature_rank = feature_rank,
    #                     feature_filter = feature_filter)
    #
    #   if (normalize) {
    #     self$normalize()
    #   }
    #   # Fetch labelled tree by featureData
    #   tree <- self$label_phylo(feature_rank = feature_rank)
    #
    #   # Subset data by correlation columns
    #   correlation_data <- na.omit(self$metaData[, .SD, .SDcols = c(sample.id, cor_columns)])
    #
    #   # Compute correlations for taxa
    #   Y <- correlation_data[, .SD, .SDcols = !c(sample.id)]
    #   cor_mat <- as.data.frame(stats::cor(x = t(as.matrix(self$countData[, correlation_data[[ sample.id ]] ])),
    #                                       y = Y,
    #                                       method = cor_method))
    #   rownames(cor_mat) <- tree$tip.label
    #   colnames(cor_mat) <- sub("CORRELATION_", "", colnames(Y))
    #
    #   # Add taxa labels of where correlation is above threshold
    #   logical <- cor_mat > cor_threshold | cor_mat < -cor_threshold
    #   logical_mat <- cor_mat[apply(logical, 1, any), ]
    #   filter_NAs <- rownames(logical_mat)[!grepl("^NA", rownames(logical_mat))]
    #
    #   # Restores omics class components
    #   private$tmp_restore()
    #
    #   # Only visuakizes taxa meeting the correlation threshold
    #   if (length(filter_NAs > 0)) {
    #     final_cor <- logical_mat[filter_NAs, ]
    #     final_tree <- ape::keep.tip(tree, tip = filter_NAs)
    #
    #     # Adding labelling layer to base tree
    #     label_offset <- length(cor_columns) * 2.5
    #     p <- ggtree(final_tree, branch.length = "none") +
    #       geom_tiplab(size = 3, offset = label_offset) +
    #       geom_treescale() +
    #       theme_tree()
    #     # Increases x-axis space based on label_offset
    #     p <- p +
    #       xlim(0, max(p$data$x) + label_offset * 2)
    #
    #     # Legend labels
    #     cor_names <- colnames(final_cor)
    #     cor_sequence <- seq_along(cor_names)
    #     column_labels <- stats::setNames(as.character(cor_sequence), cor_names)
    #
    #     # Adding heatmap to final tree
    #     return(gheatmap(p, final_cor,
    #                     offset = 0.1,
    #                     width = 1,
    #                     colnames_position = "top",
    #                     colnames_offset_y = 0.1,
    #                     custom_column_labels = cor_sequence,
    #                     hjust = 0.5,
    #                     font.size = 2.5) +
    #              scale_fill_viridis_c(option = "E",
    #                                   name = cor_method,
    #                                   na.value = "white") +
    #              labs(caption = paste(column_labels, names(column_labels),
    #                                   sep = " : ",
    #                                   collapse = " - ")))
    #   } else {
    #     return(paste0("None ", feature_rank, " are found that meet the correlation threshold of (+/-) ", cor_threshold))
    #   }
    # },
    # label_phylo = function(feature_rank) {
    #   # Create tmp tree copy
    #   tmp_tree <- self$treeData
    #
    #   # starts with empty tip labels order
    #   tip_dt <- data.table::data.table("tips" = tmp_tree$tip.label)
    #
    #   # Create lookup-table
    #   lookup_dt <- data.table::data.table("id" = self$featureData[[ "ID" ]],
    #                                       feature_rank = self$featureData[[ feature_rank ]])
    #   colnames(lookup_dt) <- c("id", feature_rank)
    #
    #   # join tables
    #   final_dt <- base::merge(tip_dt, lookup_dt, by.x="tips", by.y="id", all.x = TRUE)
    #   # Re-name tips and perform filtering if applicable.
    #   tmp_tree$tip.label <- final_dt[[ feature_rank ]]
    #
    #   return(tmp_tree)
    # },
    #' @description
    #' Automated Omics Analysis based on metadata template.
    #' For now only works with headers "RANKSTAT_" and "CORRELATION_".
    #' Samples should be as "SAMPLE-ID" upper or lower case.
    #' @param feature_ranks A character vector of features to use, default \code{c("Phylum", "Family", "Genus")}.
    #' @param feature_filter A character vector of to filter unwanted taxa, default \code{c("uncultured")}
    #' @param distance_metrics A character vector specifying what (dis)similarity metrics to use, default \code{c("unifrac")}
    #' @param dist_matrix A path to pre-computed distance matrix, expects tsv/csv/txt file from qiime2.
    #' @param alpha_div_table A path to pre-computed alpha diversity with rarefraction depth, expects tsv/csv/txt from qiime2.
    #' @param cpus Number of cores to use, only used in \link[omics]{ordination} when dist_matrix is not supplied.
    #'
    #' @return A nested list of \link[ggplot2]{ggplot} objects.
    autoFlow = function(feature_ranks = c("Phylum", "Family", "Genus"),
                        feature_filter = c("uncultured"),
                        distance_metrics = c("unifrac"),
                        dist_matrix = NULL,
                        alpha_div_table = NULL,
                        cpus = 1) {

      # Nested function
      is_empty = function(obj) {
        if (length(obj) == 0) {
          return(NULL)
        } else {
          return(obj)
        }
      }

      # Plot results as list
      plots <- list()

      # Collect columns
      metacols <- colnames(self$metaData)
      RANKSTAT_data <- self$metaData[, .SD, .SDcols = grepl("RANKSTAT_", metacols)]
      RANKSTAT_colnames <- colnames(RANKSTAT_data)
      self$metaData[, (RANKSTAT_colnames) := lapply(.SD, as.character), .SDcols = RANKSTAT_colnames]
      CORRELATION_data <- self$metaData[, .SD, .SDcols = grepl("CORRELATION_", metacols)]

      #
      #---------------------------------------------#
      # Perform standard visualizations             #
      #---------------------------------------------#
      #
      # RANKSTAT
      #
      feature_nrow <- length(feature_ranks)
      RANKSTAT_ncol <- length(RANKSTAT_data)
      #
      # Object manipulation
      #
      self$feature_subset(Kingdom == "Bacteria")
      # Standard rank stats
      plots$rankstat_plot <- self$rankstat()

      # Main loop
      if (RANKSTAT_ncol > 0) {

        # Load custom distance matrix if supplied
        if (!is.null(dist_matrix)) {
          dist_matrix <- read_tsv_matrix(filename = dist_matrix)
          dist_matrix <- dist_matrix[self$metaData[["SAMPLE-ID"]], self$metaData[["SAMPLE-ID"]]]
        }

        # Load custom rarefraction alpha diversity table if supplied
        if (!is.null(alpha_div_table)) {
          alpha_div_table <- read_rarefraction_qiime(filename = alpha_div_table)
        }

        # Initialize plot containers
        composition_plots <- matrix(list(), RANKSTAT_ncol, feature_nrow)
        correlation_plots <- list()
        Log2FC_plots <- matrix(list(), RANKSTAT_ncol, feature_nrow)
        alpha_div_plots <- list()
        metrics_nrow <- length(distance_metrics)
        pcoa_plots <- matrix(list(), RANKSTAT_ncol, metrics_nrow)
        nmds_plots <- matrix(list(), RANKSTAT_ncol, metrics_nrow)
        RDA_plots <- matrix(list(), RANKSTAT_ncol, 2)

        for (i in 1:RANKSTAT_ncol) {
          col_name <- colnames(RANKSTAT_data)[i]
          cat(paste0("Processing ... ", col_name, " \n"))

          # Alpha diversity: Shannon index
          if (inherits(alpha_div_table, "data.table")) {
            dt_final <- base::merge(alpha_div_table,
                                    self$metaData[, .SD, .SDcols = c("SAMPLE-ID", col_name)],
                                    by = "SAMPLE-ID",
                                    all.x = TRUE) %>%
              na.omit(cols = col_name)

            alpha_div_plots[[i]] <- diversity_plot(dt = dt_final,
                                                   values = "alpha_div",
                                                   col_name = col_name,
                                                   palette = fetch_palette(dt_final, col_name, "Set2"),
                                                   method = "custom")$diversity
          } else {
            alpha_div_plots[[i]] <- self$alpha_diversity(col_name = col_name, method = "shannon")$diversity
          }

          # Create RDA1 vs PC1 triplot
          RDA_plots[[i, 1]] <- self$triplot(feature_rank = "Genus",
                                            feature_filter = feature_filter,
                                            metadata.col = col_name,
                                            pairwise = FALSE,
                                            choice_dim = c("RDA1", "PC1"))$plot

          # Create PC1 vs PC2 triplot
          RDA_plots[[i, 2]] <- self$triplot(feature_rank = "Genus",
                                            feature_filter = feature_filter,
                                            metadata.col = col_name,
                                            pairwise = FALSE,
                                            choice_dim = c("PC1", "PC2"))$plot


          # Microbiome composition by all samples
          for (j in 1:feature_nrow) {
            # Creates composition long table
            res <- self$composition(feature_rank = feature_ranks[j],
                                    feature_filter = feature_filter,
                                    col_name = col_name)

            # Creates composition ggplot as list
            composition_plots[[i, j]] <- composition_plot(data = res$data,
                                                          palette = res$palette,
                                                          feature_rank = feature_ranks[j],
                                                          group_by = col_name)

            # Creates correlation ggplot as list
            if (i == 1 & length(CORRELATION_data) > 0) {
              correlation_plots[[j]] <- self$correlation(feature_rank = feature_ranks[j],
                                                         feature_filter = feature_filter,
                                                         cor_columns = colnames(CORRELATION_data))
            }


            # # Creates Log2 Fold-Change (FC) ggplot as list
            # unique_groups <- unique(na.omit(RANKSTAT_data[, .SD, .SDcols = col_name]))
            # if (nrow(unique_groups) == 2) {
            #   condition_A <- unique_groups[1, ]
            #   condition_B <- unique_groups[2, ]
            #
            #   Log2FC_plots[[i, j]] <- self$differential_feature_expression(feature_rank = feature_ranks[j],
            #                                                                sample.id = "SAMPLE-ID",
            #                                                                condition.group = col_name,
            #                                                                condition_A = condition_A,
            #                                                                condition_B = condition_B,
            #                                                                feature_filter = feature_filter)[["volcano_plot"]][[1]]
            # }


          }
          for (j in 1:metrics_nrow) {
            if (inherits(dist_matrix, "Matrix")) {
              tmp_plts <- self$ordination(distmat = dist_matrix,
                                          method = "pcoa",
                                          group_by = col_name)
            } else {
              # Creates temporary plot results for PCoA
              tmp_plts <- self$ordination(metric = distance_metrics[j],
                                          method = "pcoa",
                                          group_by = col_name,
                                          weighted = TRUE,
                                          parallel = TRUE,
                                          cpus = cpus)
            }

            pcoa_plots[[i, j]] <- patchwork::wrap_plots(tmp_plts[c("scree_plot", "anova_plot", "scores_plot")],
                                                        nrow = 1) +
              plot_layout(widths = c(rep(5, length(tmp_plts))),
                          guides = "collect")

            # Creates temporary plot results for NMDS
            if (inherits(dist_matrix, "Matrix")) {
              tmp_plts <- self$ordination(distmat = dist_matrix,
                                          method = "nmds",
                                          group_by = col_name)
            } else {
              tmp_plts <- self$ordination(metric = distance_metrics[j],
                                          method = "nmds",
                                          group_by = col_name,
                                          weighted = TRUE)
            }


            nmds_plots[[i, j]] <- patchwork::wrap_plots(tmp_plts[c("anova_plot", "scores_plot")],
                                                        nrow = 1) +
              plot_layout(widths = c(rep(5, length(tmp_plts))),
                          guides = "collect")
          }
        }
        plots$alpha_div_plots <- is_empty(alpha_div_plots)
        plots$correlation_plots <- is_empty(correlation_plots)
        plots$composition_plots <- is_empty(composition_plots)
        plots$Log2FC_plots <- is_empty(Log2FC_plots)
        plots$pcoa_plots <- is_empty(pcoa_plots)
        plots$nmds_plots <- is_empty(nmds_plots)
        plots$rda_plots <- is_empty(RDA_plots)
      }

      return(plots)
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
    eigen_80 = function(eig_explained) {
      sum_variance = 0
      counter = 1
      for (i in 1:length(eig_explained)) {
        sum_variance <- sum_variance + eig_explained[i]
        counter <- counter + 1
        if (sum_variance >= 80) break
      }

      return(counter)
    },
    subset_by_dimensions = function(model, dimensions) {
      perc_explained <- round(100*(summary(model)$cont$importance[2, dimensions]),2)
      n_dim_pairs <- dimensions[1:eigen_80(perc_explained)]
      return(perc_explained)
    },

    subset_by_species = function(model, scores_species, pc) {
      species_explained <- utils::head(base::sort(round(100*scores_species[, pc]^2, 3), decreasing = TRUE))
      scores_species_explained <- scores_species[rownames(scores_species) %in% names(species_explained),]

      result <- list(
        scores = scores_species_explained,
        explained_PC1 = species_explained
      )

      return(result)
    }
  )
)
