#' Abstract 'tools' class
#'
#' @description This is the abstract class 'tools', contains a variety of methods that are inherited and applied in the omics classes:
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @details
#' Every class is created with the \link[R6]{R6Class} method. Methods are either public or private, and only the public components are inherited by other omics classes.
#' The tools class by default uses triplet \link[data.table]{data.table} data structures for quick and efficient data manipulation and returns the object by reference, same as the R6 class.
#' The method by reference is very efficient when dealing with big data.
#' @export

tools <- R6::R6Class(
  classname = "tools",
  cloneable = FALSE,
  public = list(
    #' @field countData A path to an existing file, data.table or data.frame.
    countData = NULL,
    #' @field featureData A path to an existing file, data.table or data.frame.
    featureData = NULL,
    #' @field metaData A path to an existing file, data.table or data.frame.
    metaData = NULL,

    #' @description
    #' Wrapper function that is inherited and adapted for each omics class.
    #' To create a new object use \code{tools$new()}
    #' @param countData countData A path to an existing file, data.table or data.frame.
    #' @param featureData A path to an existing file, data.table or data.frame.
    #' @param metaData A path to an existing file, data.table or data.frame.
    #' @return A new `tools` object.
    initialize = function(countData = NA, featureData = NA, metaData = NA) {
      # Loads counts
      self$countData <- data.table::fread(countData)

      # Loads features
      self$featureData <- data.table::fread(featureData)
      self$featureData[, ID := rownames(self$featureData)]

      # Loads metadata & replaces empty values by NAs
      self$metaData <- data.table::fread(metaData)
      self$metaData <- self$metaData[, lapply(.SD, function(x) ifelse(x == "", NA, x)),
                                     .SDcols = colnames(self$metaData)]

      # Set column order
      self$countData <- self$countData[, self$metaData[["SAMPLE-ID"]], drop = FALSE]
    },
    #' @description
    #' Removes empty (zero) values by row and column.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$removeZeros()
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
    #' Remove NAs from metaData and updates the countData object fields
    #' @param sample.id The column containing sample-id that are also the countData columns.
    #' @param column The column from where NAs should be removed.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$removeNAs(sample.id = "SAMPLE-ID", column = "treatment")
    #'
    removeNAs = function(sample.id, column) {
      self$metaData <- na.omit(self$metaData, cols = column)
      self$countData <- self$countData[, self$metaData[[ sample.id ]]]
      invisible(self)
    },
    #' @description
    #' Feature subset (based on featureData), automatically applies \code{removeZeros}
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in featureData.
    #' Only rows for which all conditions evaluate to TRUE are kept.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$feature_subset(rank1 == "Streptococcus")
    #' obj$feature_subset(rank1 %in% c("Streptococcus", "uncultured"))
    feature_subset = function(...) {
      rows_to_keep <- self$featureData[, ...]
      self$featureData <- self$featureData[rows_to_keep, ]
      self$countData <- self$countData[rows_to_keep, ]
      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Sample subset (based on metaData), automatically applies \code{removeZeros}
    #' @param ... Expressions that return a logical value, and are defined in terms of the variables in metaData.
    #' Only rows for which all conditions evaluate to TRUE are kept.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$sample_subset(cycle == "t1")
    #' obj$sample_subset(cycle %in% c("t1", "t5"))
    sample_subset = function(...) {
      # set order of columns
      self$countData <- self$countData[, self$metaData[["SAMPLE-ID"]], drop = FALSE]
      # subset columns and rows
      rows_to_keep <- self$metaData[, ...]
      self$metaData <- self$metaData[rows_to_keep, ]
      self$countData <- self$countData[, rows_to_keep]
      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Agglomerates features by column, automatically applies \code{removeZeros}.
    #' @param feature_rank Column name to agglomerate.
    #' @param feature_filter Removes features by name, works on single strings or vector of strings.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$feature_glom(feature_rank = "Rank1")
    #' obj$feature_glom(feature_rank = "Genus", feature_filter = c("uncultured", "metagenome"))
    feature_glom = function(feature_rank, feature_filter = NA) {
      # creates a subset of unique feature rank, hashes combined for each unique rank
      counts <- data.table::data.table("ID" = rownames(self$countData))
      features <- data.table::copy(self$featureData)

      # set keys
      data.table::setkey(counts, ID)
      data.table::setkey(features, ID)

      # Create groups by ID
      grouped_ids <- features[, .(IDs = list(ID)), by = feature_rank]
      counts_glom <- Matrix::Matrix(0,
                                    nrow = nrow(grouped_ids),
                                    ncol = ncol(self$countData),
                                    dimnames = list(NULL, colnames(self$countData)),
                                    sparse = TRUE)

      # Populate sparse matrix by colsums of identical taxa
      for (i in 1:nrow(grouped_ids)) {
        ids <- grouped_ids$IDs[[i]]
        if (length(ids) == 1) {
          counts_glom[i, ] <- self$countData[grouped_ids$IDs[[i]],]
        } else {
          counts_glom[i, ] <- Matrix::colSums(self$countData[grouped_ids$IDs[[i]],])
        }
      }

      # Prepare final self-components
      self$featureData <- base::unique(self$featureData, by = feature_rank)
      self$countData <- counts_glom

      # Clean up featureData
      empty_strings <- self$featureData[[feature_rank]] != ""
      self$featureData <- self$featureData[empty_strings, ]
      self$countData <- self$countData[empty_strings, ]

      # Remove user-specified feature(s) filter as array
      if (is(feature_filter, "character")) {
        user_filter <- !grepl(paste(feature_filter, collapse = "|"), self$featureData[[feature_rank]])
        self$featureData <- self$featureData[user_filter, ]
        self$countData <- self$countData[user_filter, ]
      }

      self$removeZeros()
      invisible(self)
    },
    #' @description
    #' Performs transformation on countData as a Triplet sparse matrix \link[Matrix]{uniqTsparse}
    #' @param fun A function such as \code{log2}, \code{log}
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$transform(log2)
    #'
    transform = function(fun) {
      self$countData@x <- fun(self$countData@x)
      invisible(self)
    },
    #' @description
    #' Relative abundance computation by column sums.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$normalize()
    #'
    normalize = function() {
      self$countData@x <- self$countData@x / rep(Matrix::colSums(self$countData), base::diff(self$countData@p))
      invisible(self)
    },
    #---------------------------#
    # Methods for visualization #
    #---------------------------#
    #' @description
    #' Rank statistics based on featureData
    #' @details
    #' Counts the number of features identified for each column, for example in case of 16S metagenomics it would be the number of OTUs or ASVs on different taxonomy levels.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' plt <- obj$rankstat()
    #' plt
    #' @return A \link[ggplot2]{ggplot} object.
    rankstat = function() {
      # Counts number of ASVs without empty values
      values <- self$featureData[, lapply(.SD, function(x) sum(x != "")), .SDcols = !c("ID")]

      # Pivot into long table
      long_values <- data.table::melt(data = values,
                                      measure.vars = names(values),
                                      variable.name = "variable",
                                      value.name = "counts")

      # Sets order level of taxonomic ranks
      long_values[, variable := factor(variable, levels = c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain"))]


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
    #' Alpha diversity based on \link[vegan]{diversity}
    #' @param col_name The metaData column of categorical variables to create a ggplot object.
    #' @param method Diversity metric such as "shannon", "invsimpson" or "simpson"
    #' @param Brewer.palID Palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link[OmicFlow]{fetch_palette}.
    #' @param evenness A boolean wether to divide diversity by number of species, see \link[vegan]{specnumber}.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' plt <- obj$alpha_diversity(col_name = "treatment",
    #'                            method = "shannon")
    #'
    #' @return A \link[ggplot2]{ggplot} object.
    #' @seealso \link[OmicFlow]{diversity_plot}
    alpha_diversity = function(col_name, method = c("shannon", "invsimpson", "simpson"), Brewer.palID="Set2", evenness = FALSE) {

      # OUTPUT: Plot list
      plot_list <- list()

      # Save tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

      # Alpha diversity based on 'method'
      div <- data.table::data.table(diversity(x = self$countData, index=method))
      div[, (paste(col_name)) := self$metaData[, .SD, .SDcols = c(col_name)]]
      # Adjusts for evenness
      if (evenness) div$V1 <- div$V1 / log(vegan::specnumber(div$V1))

      # get colors
      colors <- fetch_palette(self$metaData, col_name, Brewer.palID)

      # Create and saves plots
      plot_list$diversity <- diversity_plot(dt = na.omit(div),
                                            values = "V1",
                                            col_name = col_name,
                                            palette = colors,
                                            method = method)

      # Restores tools class components
      private$tmp_restore()

      return(plot_list)
    },
    #' @description
    #' Visualization of compositional data.
    #' @param feature_rank A featureData column name to visualize.
    #' @param feature_filter Removes features by name, works on single strings or vector of strings.
    #' @param col_name A metaData column name to add to the compositional data.
    #' @param sample.id A metaData column specifying the sample.id, used to filter out NAs if col_name is specified.
    #' @param feature_top Integer of the top features to visualize, the max is 15, due to a limit of palettes.
    #' @param Brewer.palID Palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link[OmicFlow]{fetch_palette}.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
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
    composition = function(feature_rank, feature_filter = NA, col_name = NA, sample.id = "SAMPLE-ID", feature_top = 10, Brewer.palID = "RdYlBu") {
      # Copies object to prevent modification of tools class components
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
      if (!is.na(col_name)) {
        self$removeNAs(sample.id, col_name)
      }

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
      colnames(final_long) <- c(feature_rank, "SAMPLE-ID", "value")

      # Adds metadata columns by user input
      if (!is.na(col_name)) {
        composition_final <- base::merge(final_long,
                                         self$metaData[, .SD, .SDcols = c("SAMPLE-ID", col_name)],
                                         by = "SAMPLE-ID",
                                         all.x = TRUE)
      } else {
        composition_final <- final_long
      }


      # Factors the melted data.table by the original order of Taxa
      # Important for scale_fill_manual taxa order
      composition_final[[feature_rank]] <- factor(composition_final[[feature_rank]], levels = final_dt[[feature_rank]])

      # Restores tools class components
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
    #' @param metric A dissimilarity or similarity metric to be applied on the countData, thus far supports 'bray', 'jaccard' and 'unifrac' column name to visualize.
    #' @param method Ordination method, supports "pcoa" and "nmds".
    #' @param distmat A custom distance matrix in \link[stats]{dist} format.
    #' @param group_by A metaData column to be used as contrast for PERMANOVA or ANOSIM statistical test.
    #' @param sample.id A metaData column specifying the sample.id, used to filter out NAs from \code{group_by} column name.
    #' @param weighted Boolean, wether to compute weighted or unweighted dissimilarities.
    #' @param normalize Boolean, wether to normalize by total sample sums.
    #' @param parallel Boolean, wether to parallelize the computation of the dissimilarity matrix.
    #' @param pca.pairwise Boolean, wether to visualize different combinations of the principal components, only works with method 'pcoa'.
    #' @param pca.max.explained Integer specifying the maximum number of dissimilarity explained, used in pca.pairwise, default is 80, max number of dimensions is 15.
    #' @param pca.dim Vector with integers, specifying what dimensions to visualize in case of pca.pairwise is FALSE.
    #' @param outdir Output directory of pca.pairwise, outputs a pdf document.
    #' @param cpus Integer, number of cores to use. Default is 8 when parallelize is TRUE.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #'
    #' pcoa_plots <- obj$ordination(metric = "bray",
    #'                              method = "pcoa",
    #'                              group_by = "treatment",
    #'                              sample.id = "SAMPLE-ID",
    #'                              weighted = TRUE,
    #'                              parallel = TRUE,
    #'                              normalize = TRUE)
    #' pcoa_plots
    #'
    #' @return A list of \link[ggplot2]{ggplot} object.
    #' @seealso \link[OmicFlow]{ordination_plot}, \link[OmicFlow]{stats_plot}, \link[OmicFlow]{pairwise_anosim}, \link[OmicFlow]{pairwise_adonis}
    ordination = function(metric = c("bray", "jaccard", "unifrac"), method = c("pcoa", "nmds"), group_by, distmat = NULL, weighted = FALSE, normalize = TRUE, parallel = FALSE,
                          pca.pairwise = FALSE, pca.max.explained = 80, pca.dim = c(1,2), outdir=".", cpus = 8, sample.id = "SAMPLE-ID") {

      # Copies object to prevent modification of tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

      # Subset by missing values
      self$removeNAs(sample.id, group_by)
      if (inherits(distmat, "Matrix")) {
        distmat <- distmat[self$metaData[[sample.id]], self$metaData[[sample.id]]]
        distmat <- as.dist(distmat)
      }

      # Creates a list of plots
      plot_list <- list()

      if (parallel == TRUE) {
        # Uses available CPUs for %dopar%
        RcppParallel::setThreadOptions(numThreads = cpus)
      }

      # Normalizes counts
      if (normalize == TRUE) {
        self$normalize()
      }

      # computes distance matrix without sample rarefying
      if (is.null(distmat)) {
        # Requires rownames to contain same labels as tree
        counts <- slam::as.simple_triplet_matrix(self$countData)
        rownames(counts) <- self$featureData$ID

        if (metric == "unifrac") {
          distmat <- rbiom::beta.div(biom = counts,
                                     method = metric,
                                     weighted = weighted,
                                     tree = self$treeData)
        } else {
          distmat <- rbiom::beta.div(biom = counts,
                                     method = metric,
                                     weighted = weighted)
        }
      }

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
        "pcoa" = pairwise_adonis(distmat, groups = na.omit(self$metaData[[ group_by ]])),
        "nmds" = pairwise_anosim(distmat, groups = na.omit(self$metaData[[ group_by ]]))
      )

      # Normalization of eigenvalues
      if (method == "pcoa") {
        pcs$eig_norm <- pcs$eig %>%
          purrr::map(function(x) x / sum(pcs$eig) * 100) %>%
          unlist()

        # Collects loading scores into dataframe
        df_pcs_points <- data.table::data.table(pcs$points)
        colnames(df_pcs_points) <- base::sub("Dim", "PC", colnames(df_pcs_points))
      } else if (method == "nmds") {
        df_pcs_points <- data.table::data.table(pcs$points)
        df_pcs_points$stress <- pcs$stress
      }

      # Adds relevant data
      df_pcs_points[, groups := self$metaData[[ group_by ]] ]
      df_pcs_points[, samples := row.names(df_pcs_points) ]

      # Pairwise dimensions
      if (pca.pairwise & method == "pcoa") {
        # Finds number of dimensions that explain 80% of distances
        n_dimensions = 0
        sum_eig = 0
        for (eig in pcs$eig_norm) {
          if (sum_eig < pca.max.explained) {
            sum_eig <- sum_eig + eig
            n_dimensions <- n_dimensions + 1
          } else break
        }

        # Creates paired combinations of dimensions into a list of plots
        n_dim_pairs <- utils::combn(seq(n_dimensions), 2)
        pdf(paste0(outdir, "/pairwise_PCoA.pdf"))
        for (i in seq(ncol(n_dim_pairs))) {
          pair <- n_dim_pairs[, i]
          print(pcoa_plot(df_pcs_points, pcs, pair, metric))
        }
        dev.off()
      }

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

      # Restores tools class components
      private$tmp_restore()

      return(plot_list)
    },
    #' @description
    #' Differential feature expression
    #' @param feature_rank A featureData column name to visualize.
    #' @param feature_filter Removes features by name, works on single strings or vector of strings.
    #' @param feature_top Integer of the top features to visualize, the max is 15, due to a limit of palettes.
    #' @param sample.id A metaData column specifying the sample.id, used to filter out NAs from \code{condition.group} column name.
    #' @param paired Boolean, wether to compute paired or unpaired log2 fold change, for paired it is required to specify paired.id. Default is unpaired.
    #' @param paired.id A metaData column name containing paired ids.
    #' @param condition.group A metaData column name of where the conditions A and B are located.
    #' @param condition_A A character string or vector.
    #' @param condition_B A character string or vector.
    #' @param pvalue.threshold Integer, a P-value threshold to label and color significant features. Default is 0.05.
    #' @param foldchange.threshold Integer, a fold-change threshold to label and color significantly expressed features. Default is 0.06
    #' @param normalize Boolean, wether to normalize by total sample sums.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #'
    #' unpaired <- obj$differential_feature_expression(feature_rank = "Genus",
    #'                                                 sample.id = "SAMPLE-ID",
    #'                                                 paired = FALSE,
    #'                                                 condition.group = "treatment",
    #'                                                 condition_A = c("H"),
    #'                                                 condition_B = c("T"))
    #'
    #' paired <- obj$differential_feature_expression(feature_rank = "Genus",
    #'                                               sample.id = "SAMPLE-ID",
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
    differential_feature_expression = function(feature_rank, sample.id, paired=FALSE, paired.id,
                                               condition.group, condition_A, condition_B, pvalue.threshold=0.05, foldchange.threshold=0.06,
                                               feature_filter = NA, feature_top = NA, normalize = TRUE) {
      # Final output
      plot_list <- list()

      # Copies object to prevent modification of tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

      # Subset by missing values
      self$removeNAs(sample.id, condition.group)

      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_glom(feature_rank = feature_rank, feature_filter = feature_filter)

      # normalization if applicable
      if (normalize) {
        self$normalize()
      }

      # Check how many features to select (depended if volcano is desired)
      if (!is.na(feature_top)) {
        feature_top <- feature_top
      } else {
        feature_top <- nrow(self$featureData)
      }

      # Extract relative abundance
      rel_abun <- Matrix::rowMeans(self$countData[1:feature_top,])

      # Creates long table of relative abundance
      dt <- sparse_to_dtable(self$countData)[, (feature_rank) := self$featureData[[feature_rank]]]
      stats_dt <- base::merge(data.table::melt(dt,
                                               measure.vars = colnames(dt)[!grepl(feature_rank, colnames(dt))],
                                               variable.name = sample.id,
                                               value.name = "values"),
                              self$metaData[, .SD, .SDcols = c(sample.id, condition.group)],
                              by = sample.id)

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
      # For example A = T1, T2 and B = H1, H2
      if (paired == TRUE) {
        # sorting of metadata
        condition.labels <- data.table::setorderv(self$metaData,
                                                  cols = c(sample.id, paired.id, condition.group))[[ condition.group ]]
        # paired samples
        DFE <- paired_fold(dt = dt,
                           sample.id = sample.id,
                           paired.id = paired.id,
                           condition_A = condition_A,
                           condition_B = condition_B,
                           unique.id = unique(self$metaData[[ paired.id ]]),
                           condition_labels = condition.labels,
                           feature_rank = feature_rank,
                           cpus = cpus)
        # Save data
        plot_list$data <- DFE

      } else if (paired == FALSE) {
        # sorting of metadata
        condition.labels <- data.table::setorderv(self$metaData,
                                                  cols = c(sample.id, condition.group))[[ condition.group ]]
        # unpaired samples
        DFE <- unpaired_fold(dt = dt,
                             sample.id = sample.id,
                             condition_A = condition_A,
                             condition_B = condition_B,
                             condition_labels = condition.labels,
                             feature_rank = feature_rank,
                             cpus = cpus)

      } else {
        stop("paired can only be TRUE or FALSE, check your input.")
      }

      # Generate heatmap plot with df_diff data
      if (paired == TRUE) {
        # Adds size to paired heatmap
        add_columns <- unique(self$metaData[, .SD, .SDcols = c(sample.id, paired.id)])

        merged_data <- base::merge(
          stats_dt,
          add_columns,
          by = sample.id,
          all.x = TRUE
        )

        # Subset merged data
        subset_merged <- merged_data[, .SD, .SDcols = c(paired.id, feature_rank, "values")]
        colnames(subset_merged) <- c("SAMPLE-ID", feature_rank, "values")

        # Second merge
        final_merge <- base::merge(
          x = DFE,
          y = subset_merged,
          by = c("SAMPLE-ID", "Genus"),
          all.x = TRUE
        )

        # Check if multiple diff_ are present
        grouped_dt <- final_merge %>%
          dplyr::group_by(`SAMPLE-ID`, Genus) %>%
          dplyr::summarise(mean_values = mean(values, na.rum = TRUE),
                           diff_1 = mean(diff_1, na.rm = TRUE)) %>%
          dplyr::ungroup()


        # Generating heatmap plot based on paired boolean
        n_diff_columns <-  sum(grepl("^diff_", colnames(DFE)))

        # Generate heatmap plot with df_diff data
        heatmap_plot <- grouped_dt %>%
          ggplot(mapping = aes(x = base::get(sample.id, DFE),
                               y = base::get(feature_rank, DFE)))

        # If there is only one column uses default settings
        if (n_diff_columns == 1) {
          heatmap_plot <- heatmap_plot +
            geom_point(aes(size = mean_values, fill = diff_1), shape = 21)
        } else {
          # Adds geom_tile for number of diff_columns
          for (i in 1:n_diff_columns) {
            if (i == 1) {
              heatmap_plot <- heatmap_plot +
                geom_point(aes(size = mean_values, fill = !!sym(paste0("diff_", i))), shape = 21)
            } else {
              heatmap_plot <- heatmap_plot +
                geom_point(aes(size = mean_values, fill = !!sym(paste0("diff_", i))),
                           shape = 21,
                           position = position_nudge(x = 0.5))
            }

          }
        }
        # Finishes heatmap plot
        plot_list$tile_plot <- heatmap_plot +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
                axis.text.y = element_text(size=12),
                axis.text = element_text(size=12),
                text = element_text(size=12),
                legend.text = element_text(size=12),
                legend.title = element_text(size=14),
                strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
          scale_y_discrete(limits = rev(levels(as.factor(grouped_dt[["Genus"]])))) +
          scale_fill_gradient2(name = paste0("log2( A / B )"),
                               low = "blue",
                               mid = "white",
                               high = "red",
                               na.value = "grey80") +
          scale_size_continuous(name = "Mean Rel. Abun. (%)", labels = scales::label_number(accuracy = 0.01)) +
          labs(x = NULL,
               y = NULL)
      } else {
        #----------------------#
        # Visualization        #
        #----------------------#

        # Add relative abundance, and save data as output list
        DFE <- DFE[, "rel_abun" := rel_abun]
        plot_list$data <- DFE

        # Create & save volcano plot
        n_diff_columns <- sum(grepl("^Log2FC_", colnames(DFE)))

        plot_list$volcano_plot <- lapply(1:n_diff_columns, function(i) {
                volcano_plot(dt = DFE,
                             X = paste0("Log2FC_", i),
                             Y = paste0("pvalue_", i),
                             feature_rank = feature_rank,
                             pvalue.threshold = pvalue.threshold,
                             logfold.threshold = foldchange.threshold,
                             label_A = condition_A,
                             label_B = condition_B)
          })

      }
      # Restores tools class components
      private$tmp_restore()

      return(plot_list)
    },
    #' @description
    #' Computation and visualization of regression models
    #' Thus far it contains triplot for RDA, should be modified to perform a standard regression and then optionally also visualize RDA.
    #' @param feature_rank A featureData column name to visualize.
    #' @param feature_filter Removes features by name, works on single strings or vector of strings.
    #' @param feature_top Integer of the top features to visualize, the max is 15, due to a limit of palettes.
    #' @param metadata.col A metaData column that will be used as contrast. Multiple columns are not yet supported.
    #' @param sample.id A metaData column specifying the sample.id, used to filter out NAs from \code{metadata.col} column name.
    #' @param choice_dim A character vector to visualize which dimensions, default \code{c("RDA1", "PC1")}.
    #' @param pairwise A boolean variable, if TRUE it will perform pairwise visualization of PCA and outputs it as a pdf file.
    #' @param Brewer.palID Palette set to be applied, see \link[RColorBrewer]{brewer.pal} or \link[OmicFlow]{fetch_palette}.
    #' @param counts.scalar Adds a pseudocount to countData prior to log transformation.
    #' @return A \link[ggplot2]{ggplot} object.
    triplot = function(feature_rank, feature_filter = NA, metadata.col = NA, sample.id="SAMPLE-ID",
                       choice_dim = c("RDA1", "PC1"), pairwise = FALSE,
                       Brewer.palID = "Set2", counts.scalar = 1) {
      # Copies object to prevent modification of tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

      # Nested functions, later to be combined within tools-class
      logn <- function(otu_tab, scalar=1) {
        # log-transform, center
        # Y' = log ( A * Y + 1 ) ; where A is the 'strength' of the log transformation : 1, 10, 100, 1000, etc., default = 1
        otu_tab.log <- ( scalar * otu_tab ) + 1
        otu_tab.log <- log( otu_tab.log )
        otu_tab.sc <- scale(otu_tab.log, center = TRUE, scale = FALSE)
        return(otu_tab.sc)
      }

      eigen_80 <- function(eig_explained) {
        sum_variance = 0
        counter = 1
        for (i in 1:length(eig_explained)) {
          sum_variance <- sum_variance + eig_explained[i]
          counter <- counter + 1
          if (sum_variance >= 80) break
        }

        return(counter)
      }

      subset_by_dimensions <- function(model, dimensions) {
        perc_explained <- round(100*(summary(model)$cont$importance[2, dimensions]),2)
        n_dim_pairs <- dimensions[1:eigen_80(perc_explained)]
        return(perc_explained)
      }

      subset_by_species <- function(model, scores_species, pc) {
        species_explained <- utils::head(base::sort(round(100*scores_species[, pc]^2, 3), decreasing = TRUE))
        scores_species_explained <- scores_species[rownames(scores_species) %in% names(species_explained),]

        result <- list(
          scores = scores_species_explained,
          explained_PC1 = species_explained
        )

        return(result)
      }

      # Main pairwise code
      if (pairwise == TRUE) {
        # Creates a vector of 10 dimensions (PC1 - PC10)
        pairwise_dims <- sprintf("PC%d", seq(1:11))
        subset.result <- subset_by_dimensions(model, pairwise_dims)

        # save pdf
        pdf(paste0(outdir, "/pairwise_PCA.pdf"))
        for (i in seq(ncol(subset.result$n_dim_pairs))) {
          pair <- subset.result$n_dim_pairs[, i]
          scores_species_explained <- subset_by_species(model, scores_species, pc = pair[1])

          triplot(model, target_col, self$metaData, subset.result$var_explained, scores_species,
                  scores_species_explained, scores_sites,
                  pc1 = pair[1],
                  pc2 = pair[2])
        }
        dev.off()
      }

      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_glom(feature_rank = feature_rank, feature_filter = feature_filter)
      self$normalize()

      # Remove NAs
      if (!is.na(metadata.col)) {
        self$removeNAs(sample.id, metadata.col)
      } else stop("metadata.col is empty!")

      counts <- t(as.matrix(self$countData[, self$metaData[[ sample.id ]] ]))
      dimnames(counts)[[2]] <- self$featureData[[ feature_rank ]]

      # Subsets user specified dimensions
      pc1 <- choice_dim[1]
      pc2 <- choice_dim[2]

      # Transformation of counts and modelling to RDA
      counts.log <- logn(counts, scalar = counts.scalar)
      model <- vegan::rda(counts.log ~ get(metadata.col, self$metaData) + Condition(NULL),
                          data = self$metaData,
                          scale = FALSE,
                          na.action = na.fail,
                          subset = NULL)

      # MAIN code
      # Subset species and sites scores
      model.dims <- dim(model$CCA$u)[2] + dim(model$CA$u)[2]
      scores_species <- vegan::scores(x = model, display = "species", choices = c(1:model.dims), scaling=0)
      scores_sites <- vegan::scores(x = model, display = "sites", choices = c(1:model.dims))

      # Subset species most fitted/captured by user defined dimensions
      choice_dim.scores_species_explained <- subset_by_species(model, scores_species, pc = choice_dim[1])
      choice_dim.explained <- subset_by_dimensions(model, choice_dim)

      # Include relative abundance and significant groups in scores_sites
      rel_abun <- colSums(counts)
      Explained_species <- rownames(choice_dim.scores_species_explained$scores)
      scores_species_merged <- data.table::data.table(cbind(scores_species, rel_abun))
      scores_species_merged$taxa <- rownames(scores_species)

      # Creating color palette
      chosen_palette <- RColorBrewer::brewer.pal(length(Explained_species), Brewer.palID)
      colors <- stats::setNames(chosen_palette, Explained_species)

      # include groups for labelling and size
      scores_species_merged[, explained_species_label := ifelse(taxa %in% Explained_species, taxa, "")]
      scores_species_merged[, explained_species_size := ifelse(taxa %in% Explained_species, rel_abun, 0)]

      #Fetch groups
      mygroups <- get(metadata.col, self$metaData)

      # to be named: scores_sites
      dt <- data.table::data.table(data.frame(pc1 = scores_sites[, pc1],
                                              pc2 = scores_sites[, pc2],
                                              group = mygroups))

      # Get centroid centers for annotation
      df_mean.ord <- stats::aggregate(dt, by=list(dt$group),mean)
      colnames(df_mean.ord) <- c("Group", "x", "y")
      df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% mygroups, ]

      rslt.hull <- vegan::ordihull(scores_sites[, c(pc1, pc2)],
                                   groups = mygroups,
                                   show.group = mygroups,
                                   draw = 'none')

      # Initialize an empty data.table
      df_hull <- data.table::data.table(Group = character(), x = numeric(), y = numeric())

      # Loop through groups and bind data
      for (g in mygroups) {
        g <- as.character(g)
        x <- rslt.hull[[g]][ , 1]
        y <- rslt.hull[[g]][ , 2]
        Group <- rep(g, length(x))
        df_temp <- data.table::data.table(Group = Group, x = x, y = y)
        df_hull <- rbind(df_hull, df_temp, use.names = TRUE, fill = TRUE)
      }

      # Convert to data.table
      data.table::setDT(df_hull)

      # Restores tools class components
      private$tmp_restore()

      # Returns plot
      return(
        ggplot() +
          geom_point(data = scores_sites, aes(x = .data[[pc1]],
                                              y = .data[[pc2]]),
                     shape = 21,
                     fill = "steelblue",
                     col = "black") +
          geom_point(data = scores_species_merged, aes(x = .data[[pc1]],
                                                       y = .data[[pc2]],
                                                       size = .data[["explained_species_size"]],
                                                       col = .data[["explained_species_label"]]),
                     show.legend = TRUE,
                     stroke = ifelse(scores_species_merged$explained_species_label != "", 1.5, 0.5),
                     shape = 22) +
          geom_polygon(data=df_hull, aes(x=x, y=y, fill=Group),
                       alpha = 0.2, color = "gray40", show.legend = FALSE) +
          geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group),
                    color = "black", show.legend = FALSE) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
                axis.text.y = element_text(size=12),
                axis.text = element_text(size=12),
                text = element_text(size=12),
                legend.text = element_text(size=12),
                legend.title = element_text(size=14),
                strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
          scale_size_continuous(name = "rel. Abun.") +
          scale_color_manual(name = paste0(pc1, " explained species"),
                             values = colors) +
          labs(x = paste0(pc1, " (", choice_dim.explained[ pc1 ], "%)"),
               y = paste0(pc2, " (", choice_dim.explained[ pc2 ], "%)")) +
          guides(fill = "none", colour = guide_legend(override.aes = list(stroke = 1.5)))
      )
    },
    #' @description
    #' Correlation Analysis, default spearm and automatic filter & Visualization based on thresholds
    #' @param feature_rank A featureData column name to visualize.
    #' @param feature_filter Removes features by name, works on single strings or vector of strings.
    #' @param sample.id A metaData column specifying the sample.id, used to filter out NAs from \code{cor_columns} column name.
    #' @param cor_columns A vector of a single or multiple column names with continuous data.
    #' @param cor_method Correlation method, default spearman, see \link[stats]{cor}.
    #' @param cor_threshold Float variable, used to filter taxa that meet the cor_threshold in both positive as negative values. Default is 0.6.
    #' @param normalize Boolean, wether to normalize by total sample sums.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #'
    #' plt <- obj$correlation(feature_rank = "Genus",
    #'                        feature_filter = c("uncultured"),
    #'                        sample.id = "SAMPLE-ID",
    #'                        cor_method = "spearman",
    #'                        cor_columns = c("BMI", "weight", "flight.time"),
    #'                        cor_threshold = 0.8)
    #'
    #' @return
    #' A \link[ggplot2]{ggplot} object or \code{NULL} when no matches are found.
    correlation = function(feature_rank, feature_filter = NA, sample.id = "SAMPLE-ID",
                           cor_method = "spearman", cor_columns = c("BMI", "Weight"), cor_threshold = 0.6, normalize = TRUE) {
      # Copies object to prevent modification of tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_glom(feature_rank = feature_rank,
                        feature_filter = feature_filter)

      if (normalize) {
        self$normalize()
      }
      # Fetch labelled tree by featureData
      tree <- self$label_phylo(feature_rank = feature_rank)

      # Subset data by correlation columns
      correlation_data <- na.omit(self$metaData[, .SD, .SDcols = c(sample.id, cor_columns)])

      # Compute correlations for taxa
      Y <- correlation_data[, .SD, .SDcols = !c(sample.id)]
      cor_mat <- as.data.frame(stats::cor(x = t(as.matrix(self$countData[, correlation_data[[ sample.id ]] ])),
                                          y = Y,
                                          method = cor_method))
      rownames(cor_mat) <- tree$tip.label
      colnames(cor_mat) <- sub("CORRELATION_", "", colnames(Y))

      # Add taxa labels of where correlation is above threshold
      logical <- cor_mat > cor_threshold | cor_mat < -cor_threshold
      logical_mat <- cor_mat[apply(logical, 1, any), ]
      filter_NAs <- rownames(logical_mat)[!grepl("^NA", rownames(logical_mat))]

      # Restores tools class components
      private$tmp_restore()

      # Only visuakizes taxa meeting the correlation threshold
      if (length(filter_NAs > 0)) {
        final_cor <- logical_mat[filter_NAs, ]
        final_tree <- ape::keep.tip(tree, tip = filter_NAs)

        # Adding labelling layer to base tree
        label_offset <- length(cor_columns) * 2.5
        p <- ggtree(final_tree, branch.length = "none") +
          geom_tiplab(size = 3, offset = label_offset) +
          geom_treescale() +
          theme_tree()
        # Increases x-axis space based on label_offset
        p <- p +
          xlim(0, max(p$data$x) + label_offset * 2)

        # Legend labels
        cor_names <- colnames(final_cor)
        cor_sequence <- seq_along(cor_names)
        column_labels <- stats::setNames(as.character(cor_sequence), cor_names)

        # Adding heatmap to final tree
        return(gheatmap(p, final_cor,
                        offset = 0.1,
                        width = 1,
                        colnames_position = "top",
                        colnames_offset_y = 0.1,
                        custom_column_labels = cor_sequence,
                        hjust = 0.5,
                        font.size = 2.5) +
                 scale_fill_viridis_c(option = "E",
                                      name = cor_method,
                                      na.value = "white") +
                 labs(caption = paste(column_labels, names(column_labels),
                                      sep = " : ",
                                      collapse = " - ")))
      } else {
        return(paste0("None ", feature_rank, " are found that meet the correlation threshold of (+/-) ", cor_threshold))
      }
    },
    #' @description
    #' Relabelling phylogenetic tree by featureData
    #' @param feature_rank A featureData column name to visualize.
    #' @return A re-labelled tree as \link[ape] object
    label_phylo = function(feature_rank) {
      # Create tmp tree copy
      tmp_tree <- self$treeData

      # starts with empty tip labels order
      tip_dt <- data.table::data.table("tips" = tmp_tree$tip.label)

      # Create lookup-table
      lookup_dt <- data.table::data.table("id" = self$featureData[[ "ID" ]],
                                          feature_rank = self$featureData[[ feature_rank ]])
      colnames(lookup_dt) <- c("id", feature_rank)

      # join tables
      final_dt <- base::merge(tip_dt, lookup_dt, by.x="tips", by.y="id", all.x = TRUE)
      # Re-name tips and perform filtering if applicable.
      tmp_tree$tip.label <- final_dt[[ feature_rank ]]

      return(tmp_tree)
    },
    #' @description
    #' Automated Omics Analysis based on metadata template.
    #' For now only works with headers "RANKSTAT_" and "CORRELATION_".
    #' Samples should be as "SAMPLE-ID" upper or lower case.
    #' @param feature_ranks A character vector of features to use, default \code{c("Phylum", "Family", "Genus")}.
    #' @param feature_filter A character vector of to filter unwanted taxa, default \code{c("uncultured")}
    #' @param distance_metrics A character vector specifying what (dis)similarity metrics to use, default \code{c("unifrac")}
    #' @param dist_matrix A path to pre-computed distance matrix, expects tsv/csv/txt file from qiime2.
    #' @param alpha_div_table A path to pre-computed alpha diversity with rarefraction depth, expects tsv/csv/txt from qiime2.
    #' @param cpus Number of cores to use, only used in \link[tools]{ordination} when dist_matrix is not supplied.
    #'
    #' @return A nested list of \link[ggplot2]{ggplot} objects.
    autoFlow = function(feature_ranks = c("Phylum", "Family", "Genus"),
                        feature_filter = c("uncultured"),
                        distance_metrics = c("unifrac"),
                        dist_matrix = NA,
                        alpha_div_table = NA,
                        cpus = 1) {

      # Plot results as list
      plots <- list()

      # Collect columns
      metacols <- colnames(self$metaData)
      RANKSTAT_data <- self$metaData[, .SD, .SDcols = grepl("RANKSTAT_", metacols)]
      RANKSTAT_colnames <- colnames(RANKSTAT_data)
      self$metaData[, (RANKSTAT_colnames) := lapply(.SD, as.character), .SDcols = RANKSTAT_colnames]
      CORRELATION_data <- self$metaData[, .SD, .SDcols = grepl("CORRELATION_", metacols)]

      # Standard rank stats
      plots$rankstat_plot <- self$rankstat()
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

      # Main loop
      if (RANKSTAT_ncol > 0) {

        # Load custom distance matrix if supplied
        if (!is.na(dist_matrix)) {
          dist_matrix <- read_tsv_matrix(filename = dist_matrix)
          dist_matrix <- dist_matrix[self$metaData[["SAMPLE-ID"]], self$metaData[["SAMPLE-ID"]]]
        }

        # Load custom rarefraction alpha diversity table if supplied
        if (!is.na(alpha_div_table)) {
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
                                                   method = "custom")
          } else {
            alpha_div_plots[[i]] <- patchwork::wrap_plots(self$alpha_diversity(col_name = col_name, method = "shannon"),
                                                          nrow = 1)
          }

          # Create RDA1 vs PC1 triplot
          RDA_plots[[i, 1]] <- self$triplot(feature_rank = "Genus",
                                            feature_filter = feature_filter,
                                            metadata.col = col_name,
                                            pairwise = FALSE,
                                            choice_dim = c("RDA1", "PC1"))

          # Create PC1 vs PC2 triplot
          RDA_plots[[i, 2]] <- self$triplot(feature_rank = "Genus",
                                            feature_filter = feature_filter,
                                            metadata.col = col_name,
                                            pairwise = FALSE,
                                            choice_dim = c("PC1", "PC2"))


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
            if (i == 1) {
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

            pcoa_plots[[i, j]] <- patchwork::wrap_plots(tmp_plts,
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


            nmds_plots[[i, j]] <- patchwork::wrap_plots(tmp_plts,
                                                        nrow = 1) +
              plot_layout(widths = c(rep(5, length(tmp_plts))),
                          guides = "collect")
          }
        }
        plots$alpha_div_plots <- alpha_div_plots
        plots$correlation_plots <- correlation_plots
        plots$composition_plots <- composition_plots
        plots$Log2FC_plots <- Log2FC_plots
        plots$pcoa_plots <- pcoa_plots
        plots$nmds_plots <- nmds_plots
        plots$rda_plots <- RDA_plots
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
