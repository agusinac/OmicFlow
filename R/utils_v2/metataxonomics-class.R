metataxonomics <- R6::R6Class(
  classname = "metataxonomics",
  cloneable = FALSE,
  inherit = tools,
  public = list(
    countData = NULL,
    metaData = NULL,
    featureData = NULL,
    treeData = NULL,    
    biomData = NULL,
    initialize = function(countData = NA, metaData = NA, featureData = NA, treeData = NA, biomData = NA) {
      if (tools::file_ext(biomData) == "biom") {
        # Loads biom data
        self$biomData <- rhdf5::h5read(biomData, "/", read.attributes = TRUE)
        
        # Loads metadata
        self$metaData <- data.table::fread(metaData)
        
        # initializes count table
        private$generate_matrix()
        
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
    print = function() {
      cat("## metataxonomics-class object \n")
      if (length(self$countData) > 0) cat(paste0("## countData:\t[ ", ncol(self$countData), " Samples and ", nrow(self$countData), " Features\t] \n"))
      if (length(self$metaData) > 0) cat(paste0("## metaData:\t[ ", ncol(self$metaData), " Variables and ", nrow(self$metaData), " Samples\t] \n"))
      if (length(self$featureData) > 0) cat(paste0("## taxData:\t[ ", ncol(self$featureData)-1, " Ranks and ", nrow(self$featureData), " Taxa\t] \n"))
      if (length(self$treeData) > 0) cat(paste0("## treeData:\t[ ", length(self$treeData$tip.label), " Tips and ", self$treeData$Nnode, " Nodes\t] \n"))
    },
    reset = function() {
      self$countData = private$original_data$counts
      self$featureData = private$original_data$features
      self$metaData = private$original_data$metadata
      self$treeData = private$original_data$tree
      invisible(self)
    },
    removeZeros = function() {
      super$removeZeros()
      self$treeData <- ape::keep.tip(self$treeData, self$featureData$ID)
      invisible(self)
    }
  ),
  private = list(
    original_data = list(),
    generate_matrix = function() {
      indptr = self$biomData$sample$matrix$indptr+1
      indices = self$biomData$sample$matrix$indices+1
      data = self$biomData$sample$matrix$data
      nr = length(self$biomData$observation$ids)
      
      # Fill non-zeros among zeros
      counts = sapply(2:length(indptr), function(i) {
        x = rep(0,nr)
        seq = indptr[i-1]:(indptr[i]-1)
        x[indices[seq]] = data[seq]
        x
      })
      # save as sparse matrix
      self$countData <- data.table::data.table(counts)
      colnames(self$countData) <- self$biomData$sample$ids
    }
  )
)