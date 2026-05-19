#' Sub-class metagenomics
#'
#' @description This is a sub-class that is compatible to data obtained from either 16S rRNA marker-gene sequencing or shot-gun metagenomics sequencing.
#' It inherits all methods from the abstract class \link{omics} and only adapts the \code{initialize} function.
#' It supports BIOM format data (v2.1.0 from \url{http://biom-format.org/}) in both HDF5 and JSON format, also pre-existing data structures can be used or text files.
#' When omics data is very large, data loading becomes very expensive. It is therefore recommended to use the [`reset()`](#method-reset) method to reset your changes.
#' Every omics class creates an internal memory efficient back-up of the data, the resetting of changes is an instant process.
#' @seealso \link{omics}
#' @import R6 rhdf5 Matrix
#' @importFrom ape read.tree
#' @importFrom tools file_ext
#' @importFrom yyjsonr validate_json_file
#' @importFrom jsonlite read_json
#' @export

metagenomics <- R6::R6Class(
  classname = "metagenomics",
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
    #' Initializes the metagenomics class object with \code{metagenomics$new()}
    #' @param countData A path to an existing file or a dense/sparse \link[Matrix]{Matrix} format.
    #' @param featureData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @param metaData A path to an existing file, \link[data.table]{data.table} or data.frame.
    #' @param treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree}.
    #' @param biomData A path to an existing biom file, version 2.1.0 (http://biom-format.org/), see \link[rhdf5]{h5read}.
    #' @param feature_names A character vector to name the feature names that fit the supplied `featureData`.
    #' @return A new `metagenomics` object.
    initialize = function(countData = NULL,
                          metaData = NULL,
                          featureData = NULL,
                          treeData = NULL,
                          biomData = NULL,
                          feature_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {

      super$initialize(countData = countData,
                       featureData = featureData,
                       metaData = metaData)

      if (!is.null(biomData)) {

        if (tools::file_ext(biomData) == "biom") {
          
          #---------------------#
          ###  biomData HDF5  ###
          #---------------------#

          if (rhdf5::H5Fis_hdf5(biomData)) {

            hdf5_contents <- data.table::data.table(rhdf5::h5ls(biomData))
            hdf5_contents[, content := paste(group, name, sep = "/")]

            expected_content <- c(
              "/observation/matrix/data",
              "/observation/matrix/indptr",
              "/observation/matrix/indices",
              "/observation/ids",
              "/sample/ids")

            missing <- base::setdiff(expected_content, hdf5_contents$content)

            if (length(missing) > 0 ) {
              cli::cli_abort(
                "Expected content is missing",
                "i" = "\n{ paste(missing, collapse = ',')}"
              )
            }

            # Checks if data contains any dimensions
            list_of_dimensions <- hdf5_contents$dim[grepl(paste(expected_content, collapse="|"), hdf5_contents$content)]
            if (!all(as.numeric(list_of_dimensions) > 0)) {
              cli::cli_abort(
                "Expected content does not contain any dimensions",
                "i" = "\n{ paste(expected_content, collapse = ',')}"
              )
            }

            # Loads data in memory
            private$.biomData <- rhdf5::h5read(biomData, "/", read.attributes = TRUE)
            private$construct_hdf5_featureData()
            private$construct_hdf5_countData()

            #---------------------#
            ###  biomData JSON  ###
            #---------------------#

          } else if (yyjsonr::validate_json_file(biomData)) {
            
            private$.biomData <- jsonlite::read_json(biomData)
            private$construct_json_featureData(feature_names)
            private$construct_json_countData()

          } else {
            cli::cli_abort("{.field biomData} could not be loaded. Not a valid JSON or HDF5 format!")
          }

          # Create placeholder featureData
          if (is.null(private$.featureData)) {
            FEATURE_ID <- paste0("feature_", 1:nrow(private$.countData))
            private$.featureData <- data.table::data.table()
            private$.featureData <- private$.featureData[, (private$.feature_id) := FEATURE_ID]
            rownames(private$.countData) <- private$.featureData[[ private$.feature_id ]]
            cli::cli_alert_warning("Created placeholder {.field featureData}.")
          }
          rownames(private$.countData) <- private$.featureData[[ private$.feature_id ]]
          data.table::setcolorder(
            x = private$.featureData,
            neworder = c(private$.feature_id, base::setdiff(colnames(private$.featureData), private$.feature_id))
          )
          private$.featureData <- private$.featureData[, 
            lapply(.SD, function(x) ifelse(x == "", NA, x)),
            .SDcols = colnames(private$.featureData)]
        }
      }

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
        private$.featureData <- private$.featureData[order(match(private$.featureData[[ private$.feature_id ]], private$.treeData$tip.label))]
        private$.countData <- private$.countData[private$.featureData[[ private$.feature_id ]], ]
      }

      #-------------------#
      ###     CLEANUP   ###
      #-------------------#

      cli::cli_alert_info("Final steps .. cleaning & creating back-up")

      # Removing prefix of taxonomic features
      private$.featureData <- private$.featureData[, lapply(.SD, function(x) gsub("^[dpcofgs]_{2}", "", x)),
                                           .SDcols = colnames(private$.featureData)]
      # Rename last column names by feature_names
      n_feature_names <- length(feature_names)
      n_cols_featureData <- ncol(private$.featureData)
      colnames(private$.featureData)[n_cols_featureData:(n_cols_featureData - n_feature_names + 1)] <- base::rev(feature_names)

      # Subsetting countData by metadata
      private$sync()
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
    #' Creates a BIOM file in HDF5 format of the loaded items via ['new()'](#method-new), which is compatible to the python biom-format version 2.1, see http://biom-format.org.
    #' @param filename A character variable of either the full path of filename of the biom file (e.g. `output.biom`)
    #' @examples
    #' library("OmicFlow")
    #'
    #' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
    #' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
    #' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
    #' tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")
    #'
    #' taxa <- metagenomics$new(
    #'  metaData = metadata_file,
    #'  countData = counts_file,
    #'  featureData = features_file,
    #'  treeData = tree_file
    #' )
    #' 
    #' taxa$write_biom(filename = "output.biom")
    #' file.remove("output.biom")
    #'
    write_biom = function (filename) {

      res <- try(
        rhdf5::h5createFile(filename),
        silent = TRUE
      )
      if (!res) {
        cli::cli_abort("Can't create file {.filename {filename}}: {res}")
      }

      groups <- c(
        'observation',
        'observation/matrix',
        'observation/metadata',
        'observation/group-metadata',
        'sample',
        'sample/matrix',
        'sample/metadata',
        'sample/group-metadata'
      )

      for (group in groups)
        invisible(rhdf5::h5createGroup(filename, group))

      h5 <- try(
        rhdf5::H5Fopen(name = filename,
                       flags = 'H5F_ACC_RDWR',
                       native = TRUE),
        silent = TRUE
      )
      if (!inherits(h5, "H5IdComponent"))
        cli::cli_abort("Can't open HDF5 file {.filename {filename}}: {h5}")

      # convert countData to triplet matrix
      triplets <- Matrix::summary(private$.countData)

      #----------------------------#
      #       Add Attributes       #
      #----------------------------#
      rhdf5::h5writeAttribute(attr ="No Table ID",
                              h5obj = h5,
                              name = 'id')
      rhdf5::h5writeAttribute(attr = "OTU table",
                              h5obj = h5,
                              name = 'type')
      rhdf5::h5writeAttribute(attr = "Auto-generated biom file",
                              h5obj = h5,
                              name = 'comment')
      rhdf5::h5writeAttribute(attr = "http://biom-format.org",
                              h5obj = h5,
                              name = 'format-url')
      rhdf5::h5writeAttribute(attr = as.integer(c(2,1,0)),
                              h5obj = h5,
                              name = 'format-version')
      rhdf5::h5writeAttribute(attr = paste(Sys.Date()),
                              h5obj = h5,
                              name = 'creation-date')
      rhdf5::h5writeAttribute(attr = dim(private$.countData),
                              h5obj = h5,
                              name = 'shape')
      rhdf5::h5writeAttribute(attr = length(triplets),
                              h5obj = h5,
                              name = 'nnz')
      rhdf5::h5writeAttribute(attr = paste("OmicFlow", utils::packageVersion("OmicFlow")),
                              h5obj = h5,
                              name = 'generated-by')

      #----------------------------#
      #       Counts by row        #
      #----------------------------#
      x <- matrix(c(triplets$i - 1, triplets$j - 1, triplets$x), ncol = 3)
      x <- x[order(x[,1]),,drop=FALSE]

      counts_per_row <- base::tabulate(x[,1] + 1L, nbins = nrow(private$.countData))
      indptr <- c(0L, base::cumsum(counts_per_row))

      rhdf5::h5writeDataset(obj = base::rownames(private$.countData),
                            h5loc = h5,
                            name = 'observation/ids')
      rhdf5::h5writeDataset(obj = as.numeric(x[,3]),
                            h5loc = h5,
                            name = 'observation/matrix/data')
      rhdf5::h5writeDataset(obj = as.integer(x[,2]),
                            h5loc = h5,
                            name = 'observation/matrix/indices')
      rhdf5::h5writeDataset(obj = as.integer(indptr),
                            h5loc = h5,
                            name = 'observation/matrix/indptr')

      #----------------------------#
      #       Counts by column     #
      #----------------------------#
      x <- x[order(x[,2]),,drop=FALSE]
      counts_per_col <- base::tabulate(x[,2] + 1L, nbins = ncol(private$.countData))
      indptr <- c(0L, cumsum(counts_per_col))

      rhdf5::h5writeDataset(obj = base::colnames(private$.countData),
                            h5loc = h5,
                            name = 'sample/ids')
      rhdf5::h5writeDataset(obj = as.numeric(x[,3]),
                            h5loc = h5,
                            name = 'sample/matrix/data')
      rhdf5::h5writeDataset(obj = as.integer(x[,1]),
                            h5loc = h5,
                            name = 'sample/matrix/indices')
      rhdf5::h5writeDataset(obj = as.integer(indptr),
                            h5loc = h5,
                            name = 'sample/matrix/indptr')

      #----------------------------#
      #       Add Taxonomy         #
      #----------------------------#
      if (all(dim(private$.featureData)) > 0) {
        h5path <- 'observation/metadata/taxonomy'
        features <- as.matrix(private$.featureData[, .SD, .SDcols = !c(private$.feature_id)])
        dimnames(features) <- list(NULL, NULL)
        rhdf5::h5writeDataset(obj = features,
                              h5loc = h5,
                              name = h5path)
      }

      # Close hdf5 file connection
      rhdf5::H5Fclose(h5)
    },
    #' @description
    #' Differential feature expression (DFE) Total Sample Sum (TSS) transformed values for both paired and non-paired test.
    #' 
    #' The function performs feature agglomeration, subsetting to remove NAs in `condition.group`, finding samplepairs and finally feature scaling prior to fold-change computation. 
    #' Based on the `transform` method, fold-changes will be computed either via subtraction or division.
    #' 
    #' @param feature_rank A character or vector of characters in the `featureData` to aggregate via [`feature_merge()`](#method-feature_merge) (default: \code{"FEATURE_ID"}).
    #' @param feature_filter A character or vector of characters to remove features via regex pattern (default: \code{NULL}).
    #' @param paired A boolean value, the paired is only applicable when a `SAMPLEPAIR_ID` column exists within the `metaData`. See \link[stats]{wilcox.test} and [`samplepair_subset()`](#method-samplepair_subset).
    #' @param condition.group A character variable of an existing column name in `metaData`, wherein the conditions A and B are located.
    #' @param group_by A character variable of an existing column in `metaData` to split the table in chunks prior to fold-change computation (default: NULL).
    #' @param condition_A A character value or vector of characters.
    #' @param condition_B A character value or vector of characters.
    #' @param normalize A boolean value wether to normalize via Total Sample Sums (TSS) or not (default: \code{FALSE}).
    #' @param pvalue.threshold A numeric value used as a p-value threshold to label and color significant features (default: 0.05).
    #' @param logfold.threshold A numeric value used as a fold-change threshold to label and color significantly expressed features (default: 0.06).
    #' @param abundance.threshold A numeric value used as an abundance threshold to size the scatter dots based on their mean abundance (default: 0.01).
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
    #'  featureData = features_file
    #' )
    #'
    #' dfe <- obj$foldchange(feature_rank = "Genus",
    #'                       paired = FALSE,
    #'                       condition.group = "treatment",
    #'                       condition_A = c("healthy"),
    #'                       condition_B = c("tumor"))
    #'
    #' @returns
    #'  * `data` A long \link[data.table]{data.table} table.
    #'  * `volcano_plot` A \link[ggplot2]{ggplot} object.
    #'  * `A` A \link[data.table]{data.table} table for (each) condition A.
    #'  * `B` A \link[data.table]{data.table} table for (each) condition B.
    #'
    #' @seealso \link{volcano_plot}
    foldchange = function(
      condition.group,
      condition_A,
      condition_B,
      group_by = NULL,
      feature_rank = "FEATURE_ID",
      feature_filter = NULL,
      paired = FALSE,
      normalize = FALSE,
      pvalue.threshold = 0.05,
      logfold.threshold = 0.06,
      abundance.threshold = 0
      ) {

      ## Error handling
      #--------------------------------------------------------------------#

      if (!is.character(feature_rank) && length(feature_rank) != 1) {
        cli::cli_abort("{.val {feature_rank}} needs to be a character with a length of 1")
      } else if (!column_exists(feature_rank, private$.featureData)) {
        cli::cli_abort("The {.val {feature_rank}} column does not exist in the {.field featureData}.")
      }

      if (!is.character(condition.group) && length(condition.group) != 1) {
        cli::cli_abort("{.val {condition.group}} needs to be a character with a length of 1")
      } else if (!column_exists(condition.group, private$.metaData)) {
        cli::cli_abort("{.val {condition.group}} does not exist in the {.field metaData} or is empty.")
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

      if (!is.null(group_by)) {
        if (!is.character(group_by) && length(group_by) != 1) {
          cli::cli_abort("{.val {group_by}} needs to contain characters with length of 1.")
        } else if (!column_exists(group_by, private$.metaData)) {
          cli::cli_abort("The {.val {group_by}} column does not exist in the {.field metaData}.")
        }
      }

      ## MAIN
      #--------------------------------------------------------------------#

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

      # normalization if applicable
      if (normalize)
        self$scale(method = "tss")

      # Subset by missing values
      self$removeNAs(condition.group)

      # Subset by samplepair completion
      if (paired && !is.null(private$.samplepair_id))
        self$samplepair_subset()
      
      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_merge(feature_rank = feature_rank,
                         feature_filter = feature_filter)
      
      # Extract mean abundance
      abun <- as.matrix(Matrix::rowMeans(private$.countData))
      rownames(abun) <- private$.featureData[[ feature_rank ]]

      # Get data.table format abundances
      dt <- matrix_to_dtable(private$.countData)[, (feature_rank) := private$.featureData[[feature_rank]]]

      # Compute 2-fold expression based on (un)paired samples
      # Supports multiple inputs for A and B.
      tmp_dt <- data.table::copy(dt)

      # Apply `group_by`
      if (!is.null(group_by)) {
        chunks <- base::split(private$.metaData, by = group_by)
      } else {
        chunks <- list(all = private$.metaData)
      }
      group_names <- names(chunks)

      # subset feature labels before removing them
      feature_labels <- tmp_dt[[ feature_rank ]]
      tmp_dt <- tmp_dt[, .SD, .SDcols = !c(feature_rank)]

      # Create data.tables for results
      foldchange_dt <- data.table::data.table(feature_rank = feature_labels)
      colnames(foldchange_dt) <- feature_rank

      for (group_name in group_names) {
        condition_labels <- chunks[[group_name]][[condition.group]]

        for (i in seq_along(condition_A)) {
          # Subset by condition_A value
          dt_A <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[condition_labels %in% condition_A[i]]]
          dt_B <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[condition_labels %in% condition_B[i]]]

          # save intermediate condition tables
          output[[paste0(group_name, "_", condition_A[i])]] <- dt_A
          output[[paste0(group_name, "_", condition_B[i])]] <- dt_B

          # Convert to dense matrix
          mat_A <- as.matrix(dt_A)
          mat_B <- as.matrix(dt_B)

          # Feature means per condition
          row_means_A <- Matrix::rowMeans(mat_A)
          row_means_B <- Matrix::rowMeans(mat_B)

          # Empty vector
          result <- numeric(length(row_means_A))
          max_val <- base::max(mat_A)

          # Find zero's to prevent Inf
          both_zero <- row_means_A == 0 & row_means_B == 0
          row_means_A_zero <- row_means_A == 0 & row_means_B != 0
          row_means_B_zero <- row_means_A != 0 & row_means_B == 0
          both_non_zero <- row_means_A != 0 & row_means_B != 0

          # Compute log2 fold change
          result[both_zero] <- 0
          result[row_means_A_zero] <- row_means_A[row_means_A_zero] - log2(row_means_B[row_means_A_zero])
          result[row_means_B_zero] <- log2(row_means_A[row_means_B_zero]) - row_means_B[row_means_B_zero]
          result[both_non_zero] <- log2(row_means_A[both_non_zero]) - log2(row_means_B[both_non_zero])

          # Reverse flipped values with zero's based on max_val
          if (max_val < 1.0) {
            result[row_means_A_zero] <- result[row_means_A_zero] * -1
            result[row_means_B_zero] <- result[row_means_B_zero] * -1
          }

          # Combines to final foldchange data table
          foldchange_dt <- cbind(foldchange_dt, result)
          result_col_title <- paste0(condition_A[i], "_vs_", condition_B[i], "_in_", group_name)
          colnames(foldchange_dt)[grepl("result", colnames(foldchange_dt))] <- paste0("Log2FC_", result_col_title)

          for (k in seq_along(feature_labels)) {
            # save p-values in data.table
            suppressWarnings(
              foldchange_dt[
                k, (paste0("pvalue_", result_col_title)) := stats::wilcox.test(
                  mat_A[k, ], mat_B[k, ],
                  correct = TRUE,
                  paired = paired
                  )$p.value
                ]
            )
          }
        }
      }

      # Add abundance, and save data as output list
      dfe <- foldchange_dt[, "abun" := abun]
      output$data <- dfe

      #----------------------#
      # Visualization        #
      #----------------------#

      # Create & save volcano plot
      colnames_dfe <- colnames(dfe)
      diff_columns <- colnames_dfe[grepl("Log2FC", colnames_dfe)]
      pvalue_columns <- colnames_dfe[grepl("pvalue", colnames_dfe)]
      n_diff_columns <- length(diff_columns)

      output$volcano_plot <- lapply(1:n_diff_columns, function(i) {
        volcano_plot(
          data = dfe,
          logfold_col = diff_columns[i],
          pvalue_col = pvalue_columns[i],
          feature_rank = feature_rank,
          abundance_col = "abun",
          pvalue.threshold = pvalue.threshold,
          logfold.threshold = logfold.threshold,
          abundance.threshold = abundance.threshold,
          label_A = condition_A,
          label_B = condition_B
        ) + labs(
            subtitle = paste0(
              "Attribute: ", condition.group,
              ", test: ", ifelse(paired, "Wilcox signed rank test", "Mann-Whitney U test")
              )
          )
      })

      return(output)
    }
  ),
  private = list(
    # Private data fields
    #-------------------------#
    .countData = NULL,
    .featureData = NULL,
    .metaData = NULL,
    .treeData = NULL,
    .biomData = NULL,
    .feature_id = "FEATURE_ID",
    .sample_id = "SAMPLE_ID",
    .samplepair_id = "SAMPLEPAIR_ID",
    original_data = list(),
    construct_hdf5_featureData = function() {
      if (!is.null(private$.biomData$observation$metadata$taxonomy)) {
        private$.featureData <- data.table::data.table(t(private$.biomData$observation$metadata$taxonomy))
      }
      
      if (any(grepl(private$.feature_id, colnames(private$.metaData))) && !all(is.na(private$.metaData[[ private$.feature_id ]]))) {
        FEATURE_ID <- private$.metaData[[private$.feature_id]]
      } else if (!is.null(private$.biomData$observation$ids)) {
        FEATURE_ID <- private$.biomData$observation$ids
      } else {
        FEATURE_ID <- NULL
      }

      # Adds feature id as first column
      if (!is.null(FEATURE_ID) && !is.null(private$.featureData)) {
        private$.featureData[[ private$.feature_id ]] <- FEATURE_ID
      }

      cli::cli_alert_success("{.field featureData} is loaded.")
    },
    construct_hdf5_countData = function() {
      indptr <- as.numeric(private$.biomData$observation$matrix$indptr)

      private$.countData <- Matrix::sparseMatrix(
        i        = unlist(sapply(1:(length(indptr)-1), function (i) rep(i, diff(indptr[c(i,i+1)])))),
        j        = as.numeric(private$.biomData$observation$matrix$indices) + 1,
        x        = as.numeric(private$.biomData$observation$matrix$data),
        dims     = c(length(private$.biomData$observation$ids), length(private$.biomData$sample$ids)),
        dimnames = list(
          as.character(private$.biomData$observation$ids),
          as.character(private$.biomData$sample$ids)
        ))

      cli::cli_alert_success("{.field countData} is loaded.")
    },
    construct_json_featureData = function(feature_names) {
      # Create empty featureData
      private$.featureData <- data.table::data.table(matrix(NA_character_,
                                                     nrow = length(private$.biomData$rows),
                                                     ncol = length(c(private$.feature_id, feature_names))
                                                     ))
      setNames(private$.featureData, c(private$.feature_id, feature_names))

      # Fill first column with $id values
      private$.featureData[[ private$.feature_id ]] <- vapply(private$.biomData$rows, function(x) as.character(x$id), character(1))

      for (i in seq_along(private$.biomData$rows)) {
        taxonomy <- private$.biomData$rows[[i]]$metadata$taxonomy

        # Skip if taxonomy is missing or NULL
        if (is.null(taxonomy)) next

        # Get taxonomy index and values
        tax_indices <- seq_along(taxonomy)
        tax_values <- as.character(taxonomy)

        # Fill featureData with tax values by index
        col_positions <- tax_indices + 1
        private$.featureData[i, (col_positions) := as.list(tax_values)]
      }
      cli::cli_alert_success("{.field featureData} is loaded.")
    },
    construct_json_countData = function() {
      feature_ids <- sapply(private$.biomData$rows, function(x) unlist(x$id))
      sample_ids <- sapply(private$.biomData$columns, function(x) unlist(x$id))

      private$.countData <- Matrix::sparseMatrix(
        i        = sapply(private$.biomData$data, function(x) x[[1]]) + 1,
        j        = sapply(private$.biomData$data, function(x) x[[2]]) + 1,
        x        = sapply(private$.biomData$data, function(x) x[[3]]),
        dims     = c( length(feature_ids), length(sample_ids) ),
        dimnames = list(
          as.character(feature_ids),
          as.character(sample_ids)
        )
      )
      cli::cli_alert_success("{.field countData} is loaded.")
    }
  )
)
