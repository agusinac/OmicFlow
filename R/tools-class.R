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
      # counts
      self$countData <- data.table::fread(countData)

      # features
      self$featureData <- data.table::fread(featureData)
      self$featureData[, ID := rownames(self$featureData)]

      # metadata
      self$metaData <- data.table::fread(metaData)
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
    #' @param fun A function such as \code{function(x)}
    #' @param ... Anything following a function
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #' obj$transform(log2)
    #' obj$transform(function(x) x / sum(x))
    transform = function(fun) {
      self$countData@x <- fun(self$countData@x)
      invisible(self)
    },
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
    #' @param custom_div A custom data.frame or data.table of pre-computed diversity continuous values from qiime2 core diversity.
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
    #' plt <- obj$alpha_diversity(custom_div = shannon_df,
    #'                            col_name = "treatment")
    #' @return A \link[ggplot2]{ggplot} object.
    #' @seealso \link[OmicFlow]{diversity_plot}
    alpha_diversity = function(custom_div = NA, col_name, method = c("shannon", "invsimpson", "simpson"), Brewer.palID="Set2", evenness = FALSE) {
      # TO DO:
        # - Add hillR::hill_func_parti_pairwise(comm = counts, traits = metadata, q = 2)
        # - Find way to summarize hill results

      # OUTPUT: Plot list
      plot_list <- list(
        diversity = NULL,
        hillsPlot = NULL
      )

      # Save tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

      # Compute diversity and other metrics if custom_div is empty
      if (is.na(custom_div)) {
        # Alpha diversity based on 'method'
        div <- data.table::data.table(diversity(x = self$countData, index=method))
        div[, (paste(col_name)) := self$metaData[, .SD, .SDcols = c(col_name)]]
        # Adjusts for evenness
        if (evenness) div$V1 <- div$V1 / log(vegan::specnumber(div$V1))

        # get colors
        colors <- fetch_palette(self$metaData, col_name, Brewer.palID)

        # Create and saves plots
        plot_list$diversity <- diversity_plot(dt = div,
                                              values = "V1",
                                              col_name = col_name,
                                              palette = colors,
                                              method = method)

        # Restores tools class components
        private$tmp_restore()

        return(plot_list)

      } else {
        # Custom dataframes is used for plotting
        div <- data.table::data.table(custom_div)
        return(diversity_plot(dt = div,
                              values = div$V1,
                              col_name = div[[col_name]],
                              palette = fetch_palette(self$metaData, col_name, Brewer.palID),
                              method = method))
      }
    },
    #' @description
    #' Visualization of compositional data.
    #' @param feature_rank A featureData column name to visualize.
    #' @param feature_filter Removes features by name, works on single strings or vector of strings.
    #' @param col_name A metaData column name to add to the compositional data.
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
    composition = function(feature_rank, feature_filter = NA, col_name = NA, feature_top = 10, Brewer.palID = "RdYlBu") {
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
    #'                              weighted = TRUE,
    #'                              parallel = TRUE,
    #'                              normalize = TRUE)
    #' pcoa_plots
    #'
    #' @return A list of \link[ggplot2]{ggplot} object.
    #' @seealso \link[OmicFlow]{ordination_plot}, \link[OmicFlow]{stats_plot}, \link[OmicFlow]{pairwise_anosim}, \link[OmicFlow]{pairwise_adonis}
    ordination = function(metric = c("bray", "jaccard", "unifrac"), method = c("pcoa", "nmds"), group_by, distmat = NULL, weighted = FALSE, normalize = TRUE, parallel = FALSE,
                          pca.pairwise = FALSE, pca.max.explained = 80, pca.dim = c(1,2), outdir=".", cpus = 8) {

      # Copies object to prevent modification of tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

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
        "pcoa" = pairwise_adonis(distmat, groups = self$metaData[[ group_by ]]),
        "nmds" = pairwise_anosim(distmat, groups = self$metaData[[ group_by ]])
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

      # Creates a list of plots
      plot_list <- list(scree_plot = NULL,
                        anova_plot = NULL,
                        scores_plot = NULL)

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
    #' @param sample.id A metaData column name containing the sample ids.
    #' @param paired Boolean, wether to compute paired or unpaired log2 fold change, for paired it is required to specify paired.id. Default is unpaired.
    #' @param paired.id A metaData column name containing paired ids.
    #' @param condition.group A metaData column name of where the conditions A and B are located.
    #' @param condition_A A character string or vector.
    #' @param condition_B A character string or vector.
    #' @param pvalue.threshold Integer, a P-value threshold to label and color significant features. Default is 0.05.
    #' @param foldchange.threshold Integer, a fold-change threshold to label and color significantly expressed features. Default is 0.06
    #' @param normalize Boolean, wether to normalize by total sample sums.
    #' @param cpus Integer, number of cores to use. Default is 1.
    #' @examples
    #' obj <- tools$new(countData = "counts.csv",
    #'                  featureData = "features.txt",
    #'                  metaData = "metadata.tsv"
    #'
    #' unpaired <- obj$differential_feature_expression(feature_rank = "Genus",
    #'                                            sample.id = "SAMPLE-ID",
    #'                                            paired = FALSE,
    #'                                            condition.group = "treatment",
    #'                                            condition_A = c("H"),
    #'                                            condition_B = c("T"))
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
                                               feature_filter = NA, feature_top = NA, normalize = TRUE, cpus = 1) {
      # Final output
      plot_list <- list(
        data = NULL,
        tile_plot = NULL,
        volcano_plot = NULL
      )
      # Copies object to prevent modification of tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData,
        .treeData = self$treeData
      )

      #------#
      # Main #
      #------#

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
                              logfold.threshold = foldchange.threshold)
          })

      }
      # Restores tools class components
      private$tmp_restore()

      return(plot_list)
    },
    #' @description
    #' Computation and visualization of regression models
    regression = function() {
      # Place holder for regression, should include RDA as well
    },
    #' @description
    #' Computation and visualization of correlation models
    correlation = function(feature_rank, feature_filter = NA,
                           cor_method = "spearman", cor_columns = c("BMI", "Weight"), cor_threshold = 0.6,
                           label_offset = 10, normalize = TRUE) {
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
      correlation_data = self$metaData[, .SD, .SDcols = cor_columns]

      # Compute correlations for taxa
      cor_mat <- as.data.frame(cor(t(as.matrix(taxa$countData)), correlation_data, method = cor_method))
      rownames(cor_mat) <- tree$tip.label

      # Creating first base tree
      p <- ggtree(tree, branch.length = "none") +
        geom_tiplab(size = 3,
                    offset = label_offset) +
        geom_tippoint() +
        geom_treescale() +
        theme_tree()

      # Add taxa labels of where correlation is above threshold
      rows_to_keep <- apply(cor_mat > cor_threshold | cor_mat < -cor_threshold, 1, any)
      feature_tippoints <- rownames(cor_mat[rows_to_keep ,])
      tip_labels <- p$data[p$data$label %in% feature_tippoints, ]

      # Adding labelling layer to base tree
      p1 <- p +
        geom_tippoint(data = tip_labels,
                      mapping = aes(x = x,
                                    y = y,
                                    label = label),
                      color = "red") +
        scale_color_manual(labels = c("black" = "weak",
                                      "red" = "strong"),
                           name = "Correlation strength") +
        xlim(0, max(p$data$x) + label_offset * 4)

      # Restores tools class components
      private$tmp_restore()

      # Adding heatmap to final tree
      return(
        gheatmap(p1, cor_mat,
                 offset = 0.1,
                 width = 0.3,
                 colnames_position = "top",
                 hjust = 0.5,
                 font.size = 2.5) +
          scale_fill_viridis_c(option = "E",
                               name = cor_method,
                               na.value = "white")
      )
    },
    #' @description
    #' Relabelling phylogenetic tree by featureData
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
    #' @param feature_ranks A character vector of features to use.
    #' @param distance_metrics A character vector of dissimilarity metrics to use.
    #' @param output String variable of the out folder.
    #' @param shannon_table A path to pre-computed alpha diversity file
    #' @param distance_matrix A path to pre-computed distance matrix
    #'
    #' @return A nested list of \link[ggplot2]{ggplot} objects.
    autoFlow = function(feature_ranks = c("Phylum", "Family", "Genus", "Species"), distance_metrics = c("unifrac","bray"), output = NA, shannon_table, distance_matrix) {
      # Plot results as list
      plots <- list(
        rankstat_plot = NULL,
        shannon_plots = NULL,
        pcoa_plots = NULL,
        nmds_plots = NULL,
        composition_plots = NULL,
        correlation_heatmap_plt = NULL,
        heatmap_plots = NULL,
        RDA_plots = NULLs
      )

      # Collect columns
      metacols <- colnames(self$metaData)

      RANKSTAT_data <- self$metaData[, .SD, .SDcols = grepl("RANKSTAT_", metacols)]
      CORRELATION_data <- self$metaData[, .SD, .SDcols = grepl("CORRELATION_", metacols)]
      PAIREDGROUPBY_data <- self$metaData[, .SD, .SDcols = grepl("PAIREDGROUPBY_", metacols)]

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
      self$feature_subset(Domain == "Bacteria")
      self$transform(function(x) x / sum(x))

      # Main loop
      if (RANKSTAT_ncol > 0) {

        composition_plots <- matrix(list(), RANKSTAT_ncol, feature_nrow)
        shannon_plots <- list()
        metrics_nrow <- length(metrics)
        pcoa_plots <- matrix(list(), RANKSTAT_ncol, nrow)
        nmds_plots <- matrix(list(), RANKSTAT_ncol, nrow)

        for (i in 1:RANKSTAT_ncol) {
          col_name <- colnames(RANKSTAT_data)[i]

          # Alpha diversity: Shannon index
          shannon_plots[[i]] <- self$shannon(df_shannon = data.table::data.table(shannon_table),
                                             col_name = col_name)

          # Microbiome composition by all samples
          for (j in 1:feature_nrow) {
            # Creates composition long table
            res <- self$composition(feature_rank = feature_ranks[j],
                                    feature_filter = c("uncultured"))

            # Creates composition ggplot as list
            composition_plots[[i, j]] <- composition_plot(data = res$data,
                                                          palette = res$palette,
                                                          feature_rank = feature_ranks[j])
          }
          for (j in 1:metrics_nrow) {
            pcoa_plots[[i, j]] <- patchwork::wrap_plots(self$ordination(metric = metrics[j],
                                                                        method = "pcoa",
                                                                        weighted = TRUE),
                                                        nrow = 1) +
              plot_layout(widths = c(5, 5, 5),
                          guides = "collect")

            nmds_plots[[i, j]] <- patchwork::wrap_plots(self$ordination(metric = metrics[j],
                                                                        method = "nmds",
                                                                        weighted = TRUE),
                                                        nrow = 1) +
              plot_layout(widths = c(5, 5, 5),
                          guides = "collect")
          }
        }
        plots$shannon_plots <- shannon_plots
        plots$composition_plots <- composition_plots
        plots$pcoa_plots <- pcoa_plots
        plots$nmds_plots <- nmds_plots
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
