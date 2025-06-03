#' Sub-class metataxonomics
#'
#' @description This is a sub-class for 16S metagenomics data, called metataxonomics.
#' It inherits all methods from the abstract class \link[OmicFlow]{tools} and only adapts the \code{initialize} function.
#'
#' @export

metataxonomics <- R6::R6Class(
  classname = "metataxonomics",
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
    #' Initializes the metataxonomics class object with \code{metataxonomics$new()}
    #' @param countData countData A path to an existing file or sparseMatrix.
    #' @param featureData A path to an existing file, data.table or data.frame.
    #' @param metaData A path to an existing file, data.table or data.frame.
    #' @param treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree}.
    #' @param biomData A path to an existing biom file, version 2.1.0, see \link[rhdf5]{h5read}.
    #' @examples
    #' taxa <- metataxonomics$new(metaData = "metadata.tsv",
    #'                            biomData = "biom_with_taxonomy.biom",
    #'                            treeData = "rooted_tree.newick")
    #'
    #' @return A new `metataxonomics` object.
    initialize = function(countData = NA, metaData = NA, featureData = NA, treeData = NA, biomData = NA,
                          feature_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      if (!is.na(biomData) & tools::file_ext(biomData) == "biom") {
        # Loads biom data
        self$biomData <- rbiom::read.biom(biomData)

        # Loads metadata & replaces empty values by NAs
        self$metaData <- data.table::fread(metaData, header=TRUE)
        self$metaData <- self$metaData[, lapply(.SD, function(x) ifelse(x == "", NA, x)),
                                       .SDcols = colnames(self$metaData)]

        # Converts slam::triple_sparse_matrix to Matrix::sparseMatrix
        self$countData <- Matrix::sparseMatrix(i = self$biomData$counts$i,
                                               j = self$biomData$counts$j,
                                               x = self$biomData$counts$v,
                                               dimnames = self$biomData$counts$dimnames)

        # Match and order row names for countData and metaData
        # Now assuming countData has the right order and mistake is made in metadata
        self$countData <- self$countData[, colnames(self$countData) %in% self$metaData[["SAMPLE-ID"]], drop = FALSE]
        self$metaData <- self$metaData[self$metaData[["SAMPLE-ID"]] %in% colnames(self$countData), ]
        self$countData <- self$countData[, self$metaData[["SAMPLE-ID"]], drop = FALSE]

        # initializes taxonomy table
        if (!is.null(self$biomData$taxonomy)) {
          self$featureData <- data.table::data.table(self$biomData$taxonomy)
          colnames(self$featureData) <- sub(".otu", "ID", colnames(self$biomData$taxonomy))
        }

      } else {
        super$initialize(countData = countData,
                         featureData = featureData,
                         metaData = metaData)
      }

      if (!is.null(self$featureData)) {
        self$featureData <- self$featureData[, lapply(.SD, function(x) gsub("^[dpcofgs]_{2}", "", x)),
                                             .SDcols = colnames(self$featureData)]
        # Rename last column names by feature_names
        n_feature_names <- length(feature_names)
        n_cols_featureData <- ncol(self$featureData)
        colnames(self$featureData)[n_cols_featureData:(n_cols_featureData - n_feature_names + 1)] <- base::rev(feature_names)
      }

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
    #' Displays parameters of the `metataxonomics` object via stdout.
    #' @examples
    #' taxa <- metataxonomics$new(metaData = "metadata.tsv",
    #'                            biomData = "biom_with_taxonomy.biom",
    #'                            treeData = "rooted_tree.newick")
    #'
    #' # method 1 to call print function
    #' taxa
    #'
    #' # method 2 to call print function
    #' taxa$print()
    #'
    print = function() {
      cat("## metataxonomics-class object \n")
      if (length(self$countData) > 0) cat(paste0("## countData:\t[ ", ncol(self$countData), " Samples and ", nrow(self$countData), " Features\t] \n"))
      if (length(self$metaData) > 0) cat(paste0("## metaData:\t[ ", ncol(self$metaData), " Variables and ", nrow(self$metaData), " Samples\t] \n"))
      if (length(self$featureData) > 0) cat(paste0("## taxData:\t[ ", ncol(self$featureData)-1, " Ranks and ", nrow(self$featureData), " Taxa\t] \n"))
      if (length(self$treeData) > 0) cat(paste0("## treeData:\t[ ", length(self$treeData$tip.label), " Tips and ", self$treeData$Nnode, " Nodes\t] \n"))
    },
    #' @description
    #' Upon creation of a new `metataxonomics` object a small backup of the original data is created.
    #' Since modification of the object is done by reference and duplicates are not made, it is possible to `reset` changes to the class.
    #' The methods from the abstract class `tools` also contain a private method to prevent any changes to the original object. Such cases are ordination, alpha_diversity, differential_feature_expression.
    #' @examples
    #' taxa <- metataxonomics$new(metaData = "metadata.tsv",
    #'                            biomData = "biom_with_taxonomy.biom",
    #'                            treeData = "rooted_tree.newick")
    #'
    #' # Performs modifications
    #' taxa$transform(log2)
    #'
    #' # resets
    #' taxa$reset()
    #'
    #' # An inbuilt reset function prevents unwanted modification to the taxa object.
    #' taxa$rankstat()
    #'
    reset = function() {
      self$countData = private$original_data$counts
      self$featureData = private$original_data$features
      self$metaData = private$original_data$metadata
      self$treeData = private$original_data$tree
      invisible(self)
    },
    #' @description
    #' Removes empty (zero) values by row, column and tips.
    #' @examples
    #' taxa <- metataxonomics$new(metaData = "metadata.tsv",
    #'                            biomData = "biom_with_taxonomy.biom",
    #'                            treeData = "rooted_tree.newick")
    #'
    #' # Sample subset induces empty features
    #' taxa$sample_subset(cycle == "t1")
    #'
    #' # Remove empty features from countData and treeData
    #' taxa$removeZeros()
    removeZeros = function() {
      super$removeZeros()
      if (!is.null(self$treeData)) self$treeData <- ape::keep.tip(self$treeData, self$featureData$ID)
      invisible(self)
    },
    #' @description
    #' Writes \code{metataxonomics$new()} object as biom file, in a biom-format version 2.1 compatible format.
    #' @param file Filename of the output biom file
    #' @examples
    #' taxa <- metataxonomics$new(metaData = "metadata.tsv",
    #'                            biomData = "biom_with_taxonomy.biom",
    #'                            treeData = "rooted_tree.newick")
    #'
    #' taxa$write_biom(file = "new_output.biom")
    #'
    write_biom = function(file) {
      # Create empty biom file
      rhdf5::h5createFile(file)
      # Create groups
      invisible(rhdf5::h5createGroup(file = file, group = '/observation'))
      invisible(rhdf5::h5createGroup(file = file, group = '/observation/matrix'))
      invisible(rhdf5::h5createGroup(file = file, group = '/observation/metadata'))
      invisible(rhdf5::h5createGroup(file = file, group = '/observation/group-metadata'))
      invisible(rhdf5::h5createGroup(file = file, group = '/sample'))
      invisible(rhdf5::h5createGroup(file = file, group = '/sample/matrix'))
      invisible(rhdf5::h5createGroup(file = file, group = '/sample/metadata'))
      invisible(rhdf5::h5createGroup(file = file, group = '/sample/group-metadata'))

      # Read file
      h5 <- rhdf5::H5Fopen(file)

      # Add Attributes
      rhdf5::h5writeAttribute(attr = paste("No Table ID"),
                              h5obj = h5,
                              name = 'id')
      rhdf5::h5writeAttribute(attr = paste("Table"),
                              h5obj = h5,
                              name = 'type')
      rhdf5::h5writeAttribute(attr = "http://biom-format.org",
                              h5obj = h5,
                              name = 'format-url')
      rhdf5::h5writeAttribute(attr = as.integer(c(2,1,0)),
                              h5obj = h5,
                              name = 'format-version',
                              encoding = 3)
      rhdf5::h5writeAttribute(attr = paste(Sys.Date()),
                              h5obj = h5,
                              name = 'creation-date')
      rhdf5::h5writeAttribute(attr = dim(self$countData),
                              h5obj = h5,
                              name = 'shape',
                              encoding = 2)
      rhdf5::h5writeAttribute(attr = max(self$countData@p),
                              h5obj = h5,
                              name = 'nnz')
      rhdf5::h5writeAttribute(attr = paste("OmicFlow", utils::packageVersion("OmicFlow")),
                              h5obj = h5,
                              name = 'generated-by')

      # Read counts by taxa
      x <- matrix(c(self$countData@i - 1, self$countData@p - 1, self$countData@x), byrow=FALSE, ncol=3)

      x <- x[order(x[,1]),,drop=FALSE]
      indptr <- cumsum(unname(table(factor(x[,1]+1, 0:nrow(self$countData)))))

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

      # Read counts by sample
      x <- x[order(x[,2]),,drop=FALSE]
      indptr <- base::cumsum(base::unname(base::table(base::factor(x[,2]+1, 0:ncol(self$countData)))))

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

      # Add Taxonomy
      if (ncol(self$featureData) > 1) {
        h5path <- 'observation/metadata/taxonomy'
        features <- t(as.matrix(self$featureData[, .SD, .SDcols = !c("ID")]))
        dimnames(features) <- list(NULL, NULL)
        rhdf5::h5writeDataset(obj = features,
                              h5loc = h5,
                              name = h5path)
      }
      # Close biom file connection
      rhdf5::H5Fflush(h5)
      rhdf5::H5Fclose(h5)
    }
  ),
  private = list(
    original_data = list()
    )
)
