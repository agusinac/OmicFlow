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

    #' @field biomData A path to an existing biom file or hdf5 file, see \link[rhdf5]{h5read}.
    biomData = NULL,

    #' @description
    #' Initializes the metagenomics class object with \code{metagenomics$new()}
    #' @param countData countData A path to an existing file or sparseMatrix.
    #' @param featureData A path to an existing file, data.table or data.frame.
    #' @param metaData A path to an existing file, data.table or data.frame.
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
            self$biomData <- rhdf5::h5read(biomData, "/", read.attributes = TRUE)
            private$construct_hdf5_featureData()
            private$construct_hdf5_countData()

            #---------------------#
            ###  biomData JSON  ###
            #---------------------#

          } else if (yyjsonr::validate_json_file(biomData)) {
            
            self$biomData <- jsonlite::read_json(biomData)
            private$construct_json_featureData(feature_names)
            private$construct_json_countData()

          } else {
            cli::cli_abort("biomData could not be loaded. Not a valid JSON or HDF5 format!")
          }
        }
      }

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

      #-------------------#
      ###     CLEANUP   ###
      #-------------------#

      cli::cli_alert_info("Final steps .. cleaning & creating back-up")

      # Removing prefix of taxonomic features
      self$featureData <- self$featureData[, lapply(.SD, function(x) gsub("^[dpcofgs]_{2}", "", x)),
                                           .SDcols = colnames(self$featureData)]
      # Rename last column names by feature_names
      n_feature_names <- length(feature_names)
      n_cols_featureData <- ncol(self$featureData)
      colnames(self$featureData)[n_cols_featureData:(n_cols_featureData - n_feature_names + 1)] <- base::rev(feature_names)

      # Subsetting countData by metadata
      self$countData <- self$countData[, self$metaData[[ self$.sample_id ]], drop = FALSE]

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
    #' Displays parameters of the metagenomics object via stdout.
    #' @examples
    #' taxa_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' taxa <- readRDS(taxa_path)
    #'
    #' # method 1 to call print function
    #' taxa
    #'
    #' # method 2 to call print function
    #' taxa$print()
    #'
    #' @return object in place
    print = function() {
      cat("## metagenomics-class object \n")
      if (length(self$countData) > 0) cat(paste0("## countData:\t[ ", ncol(self$countData), " Samples and ", nrow(self$countData), " Features\t] \n"))
      if (length(self$metaData) > 0) cat(paste0("## metaData:\t[ ", ncol(self$metaData), " Variables and ", nrow(self$metaData), " Samples\t] \n"))
      if (length(self$featureData) > 0) cat(paste0("## taxData:\t[ ", ncol(self$featureData)-1, " Ranks and ", nrow(self$featureData), " Taxa\t] \n"))
      if (length(self$treeData) > 0) cat(paste0("## treeData:\t[ ", length(self$treeData$tip.label), " Tips and ", self$treeData$Nnode, " Nodes\t] \n"))
    },
    #' @description
    #' Upon creation of a new `metagenomics` object a small backup of the original data is created.
    #' Since modification of the object is done by reference and duplicates are not made, it is possible to `reset` changes to the class.
    #' The methods from the abstract class \link{omics} also contains a private method to prevent any changes to the original object when using methods such as \code{ordination} \code{alpha_diversity} or \code{$DFE}.
    #' @examples
    #' library(ggplot2)
    #'
    #' taxa_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' taxa <- readRDS(taxa_path)
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
    #' taxa_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' taxa <- readRDS(taxa_path)
    #' 
    #' # Sample subset induces empty features
    #' taxa$sample_subset(treatment == "tumor")
    #'
    #' # Remove empty features from countData and treeData
    #' taxa$removeZeros()
    #' 
    #' @return object in place
    removeZeros = function() {
      super$removeZeros()
      if (!is.null(self$treeData)) self$treeData <- ape::keep.tip(self$treeData, self$featureData$FEATURE_ID)
      invisible(self)
    },
    #' @description
    #' Creates a BIOM file in HDF5 format of the loaded items via ['new()'](#method-new), which is compatible to the python biom-format version 2.1, see http://biom-format.org.
    #' @param filename A character variable of either the full path of filename of the biom file (e.g. `output.biom`)
    #' @examples
    #' taxa_path <- system.file("extdata", "mock_taxa.rds", package = "OmicFlow", mustWork = TRUE)
    #' taxa <- readRDS(taxa_path)
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
      triplets <- Matrix::summary(self$countData)

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
      rhdf5::h5writeAttribute(attr = dim(self$countData),
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

      counts_per_row <- base::tabulate(x[,1] + 1L, nbins = nrow(self$countData))
      indptr <- c(0L, base::cumsum(counts_per_row))

      rhdf5::h5writeDataset(obj = base::rownames(self$countData),
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
      counts_per_col <- base::tabulate(x[,2] + 1L, nbins = ncol(self$countData))
      indptr <- c(0L, cumsum(counts_per_col))

      rhdf5::h5writeDataset(obj = base::colnames(self$countData),
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
      if (all(dim(self$featureData)) > 0) {
        h5path <- 'observation/metadata/taxonomy'
        features <- as.matrix(self$featureData[, .SD, .SDcols = !c("FEATURE_ID")])
        dimnames(features) <- list(NULL, NULL)
        rhdf5::h5writeDataset(obj = features,
                              h5loc = h5,
                              name = h5path)
      }

      # Close hdf5 file connection
      rhdf5::H5Fclose(h5)
    }
  ),
  private = list(
    original_data = list(),
    construct_hdf5_featureData = function() {
      self$featureData <- data.table::data.table(t(self$biomData$observation$metadata$taxonomy))

      if (any(grepl(self$.feature_id, colnames(self$metaData))) && !all(is.na(self$metaData[[ self$.feature_id ]]))) {
        FEATURE_ID <- self$metaData[[self$.feature_id]]
      } else {
        FEATURE_ID <- self$biomData$observation$ids
      }

      # Adds feature id as first column
      self$featureData[[ self$.feature_id ]] <- FEATURE_ID
      data.table::setcolorder(x = self$featureData,
                              neworder = c(self$.feature_id, base::setdiff(colnames(self$featureData), self$.feature_id))
                              )
      colnames(self$featureData) <- gsub("\\s+", "_", colnames(self$featureData))
      self$featureData <- self$featureData[, lapply(.SD, function(x) ifelse(x == "", NA, x)),
                                           .SDcols = colnames(self$featureData)]

      cli::cli_alert_success("featureData is loaded.")
    },
    construct_hdf5_countData = function() {
      indptr <- as.numeric(self$biomData$observation$matrix$indptr)

      self$countData <- Matrix::sparseMatrix(
        i        = unlist(sapply(1:(length(indptr)-1), function (i) rep(i, diff(indptr[c(i,i+1)])))),
        j        = as.numeric(self$biomData$observation$matrix$indices) + 1,
        x        = as.numeric(self$biomData$observation$matrix$data),
        dims     = c(length(self$biomData$observation$ids), length(self$biomData$sample$ids)),
        dimnames = list(
          as.character(self$biomData$observation$ids),
          as.character(self$biomData$sample$ids)
        ))

      rownames(self$countData) <- self$featureData$FEATURE_ID

      cli::cli_alert_success("countData is loaded.")
    },
    construct_json_featureData = function(feature_names) {
      # Create empty featureData
      self$featureData <- data.table::data.table(matrix(NA_character_,
                                                        nrow = length(self$biomData$rows),
                                                        ncol = length(c(self$.feature_id, feature_names))))
      setNames(self$featureData, c(self$.feature_id, feature_names))

      # Fill first column with $id values
      self$featureData[["FEATURE_ID"]] <- vapply(self$biomData$rows, function(x) as.character(x$id), character(1))

      for (i in seq_along(self$biomData$rows)) {
        taxonomy <- self$biomData$rows[[i]]$metadata$taxonomy

        # Skip if taxonomy is missing or NULL
        if (is.null(taxonomy)) next

        # Get taxonomy index and values
        tax_indices <- seq_along(taxonomy)
        tax_values <- as.character(taxonomy)

        # Fill featureData with tax values by index
        col_positions <- tax_indices + 1
        self$featureData[i, (col_positions) := as.list(tax_values)]
      }

      if (any(grepl(self$.feature_id, colnames(self$metaData))) && !all(is.na(self$metaData[[self$.feature_id]]))) {
        FEATURE_ID <- self$metaData[[self$.feature_id]]
      }

      # Adds feature id as first column
      self$featureData[, (self$.feature_id) := FEATURE_ID]
      data.table::setcolorder(x = self$featureData,
                              neworder = c(self$.feature_id, base::setdiff(colnames(self$featureData), self$.feature_id))
      )
      colnames(self$featureData) <- gsub("\\s+", "_", colnames(self$featureData))

      cli::cli_alert_success("featureData is loaded.")
    },
    construct_json_countData = function() {
      feature_ids <- sapply(self$biomData$rows, function(x) unlist(x$id))
      sample_ids <- sapply(self$biomData$columns, function(x) unlist(x$id))

      self$countData <- Matrix::sparseMatrix(
        i        = sapply(self$biomData$data, function(x) x[[1]]) + 1,
        j        = sapply(self$biomData$data, function(x) x[[2]]) + 1,
        x        = sapply(self$biomData$data, function(x) x[[3]]),
        dims     = c( length(feature_ids), length(sample_ids) ),
        dimnames = list(
          as.character(feature_ids),
          as.character(sample_ids)
        )
      )

      rownames(self$countData) <- self$featureData$FEATURE_ID

      cli::cli_alert_success("countData is loaded.")
    }
  )
)
