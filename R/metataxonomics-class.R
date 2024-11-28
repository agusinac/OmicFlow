#' Sub-class metataxonomics
#'
#' @description This is a sub-class for 16S metagenomics data, called metataxonomics.
#' It inherits all methods from the abstract class \link[OmicFlow]{tools} and only adapts the \code{initialize} function.
#'
#' @export

metataxonomics <- R6::R6Class(
  classname = "metataxonomics",
  cloneable = FALSE,
  inherit = tools,
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
    #' @param countData countData A path to an existing file, data.table or data.frame.
    #' @param featureData A path to an existing file, data.table or data.frame.
    #' @param metaData A path to an existing file, data.table or data.frame.
    #' @param treeData A path to an existing newick file or class "phylo", see \link[ape]{read.tree}.
    #' @param biomData A path to an existing biom file or hdf5 file, see \link[rhdf5]{h5read}.
    #' @examples
    #' taxa <- metataxonomics$new(metaData = "metadata.tsv",
    #'                            biomData = "biom_with_taxonomy.biom",
    #'                            treeData = "rooted_tree.newick")
    #'
    #' @return A new `metataxonomics` object.
    initialize = function(countData = NA, metaData = NA, featureData = NA, treeData = NA, biomData = NA) {
      if (tools::file_ext(biomData) == "biom") {
        # Loads biom data
        self$biomData <- rhdf5::h5read(biomData, "/", read.attributes = TRUE)

        # Loads metadata
        self$metaData <- data.table::fread(metaData)

        # initializes count matrix
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

        # Set column order
        self$countData <- self$countData[, self$metaData[["SAMPLE-ID"]], drop = FALSE]

        # initializes taxonomy table
        self$featureData <- data.table::data.table(t(self$biomData$observation$metadata$taxonomy))
        self$featureData <- self$featureData[, ID := self$biomData$observation$ids]

      } else {
        super$initialize()
      }

      colnames(self$featureData)[!grepl("ID", colnames(self$featureData))] <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      self$featureData <- self$featureData[, lapply(.SD, function(x) gsub("^[dpcofgs]_{2}", "", x)),
                                           .SDcols = colnames(self$featureData)]

      self$treeData <- ape::read.tree(treeData)

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
      self$treeData <- ape::keep.tip(self$treeData, self$featureData$ID)
      invisible(self)
    }
  ),
  private = list(
    original_data = list()
    )
)
