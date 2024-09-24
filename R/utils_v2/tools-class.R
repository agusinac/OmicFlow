tools <- R6::R6Class(
  classname = "tools",
  cloneable = FALSE,
  public = list(
    countData = NULL,
    featureData = NULL,
    metaData = NULL,
    initialize = function(countData = NA, featureData = NA, metaData = NA) {
      # counts
      self$countData <- data.table::fread(countData)
      
      # features
      self$featureData <- data.table::fread(featureData)
      self$featureData[, ID := rownames(self$featureData)]
      
      # metadata
      self$metaData <- data.table::fread(metaData)
    },
    #----------------------------#
    # Methods for data wrangling #
    #----------------------------#
    removeZeros = function() {
      # Remove empty samples (columns)
      keep_cols <- self$countData[, lapply(.SD, sum) > 0,
                                  .SDcols = colnames(self$countData)] 
      
      # Remove empty species (rows)
      keep_rows <- self$countData[, self$countData[, sum(.SD) > 0, 
                                                   .SDcols = colnames(self$countData), 
                                                   by = .I]$V1]
      # Creates new countData instance
      self$countData <- self$countData[keep_rows, .SD, .SDcols = keep_cols]
      self$featureData <- self$featureData[keep_rows]
      invisible(self)
    },
    feature_subset = function(...) {
      rows_to_keep <- self$featureData[, ...]
      self$featureData <- self$featureData[rows_to_keep, ]
      self$countData <- self$countData[rows_to_keep, ]
      self$removeZeros()
      invisible(self)
    },
    sample_subset = function(...) {
      # set order of columns
      data.table::setcolorder(self$countData, self$metaData$`SAMPLE-ID`)
      # subset columns and rows
      rows_to_keep <- self$metaData[, ...]
      self$metaData <- self$metaData[rows_to_keep, ]
      self$countData <- self$countData[, .SD, .SDcols = rows_to_keep]
      self$removeZeros()
      invisible(self)
    },
    feature_glom = function(feature_rank, feature_filter = NA) {
      # creates a subset of unique feature rank, hashes combined for each unique rank
      id_list <- self$featureData[, ID]
      counts <- self$countData[, ID := id_list]
      features <- self$featureData
      
      # set keys
      data.table::setkey(counts, ID)
      data.table::setkey(features, ID)
      
      grouped_ids <- features[, .(IDs = list(ID)), by = feature_rank]
      list_counts <- lapply(grouped_ids$ID, function(id) {
        counts[id, colSums(.SD), .SDcols = !c("ID"), on = "ID"]
      })
      self$countData <- data.table::data.table(do.call("rbind", list_counts))
      
      # subset feature data
      self$featureData <- base::unique(self$featureData, by = feature_rank)
      
      # Remove empty strings
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
    transform = function(fun, ...) {
      tmp_trans <- apply(as(as.matrix(self$countData), "TsparseMatrix"), 2, fun, ...)
      tmp_trans[!is.finite(tmp_trans)] <- 0
      self$countData <- data.table::data.table(tmp_trans)
      invisible(self)
    },
    #---------------------------#
    # Methods for visualization #
    #---------------------------#
    rankstat = function() {
      # Counts number of ASVs without empty values
      values <- self$featureData[, lapply(.SD, function(x) sum(x != "")), .SDcols = !c("ID")]
      
      # Pivot into long table
      long_values <- data.table::melt(data = values, 
                                      measure.vars = names(values), 
                                      variable.name = "variable", 
                                      value.name = "counts")
      
      # Returns rankstat plot
      return(long_values %>% 
               ggplot(mapping = aes(x = base::rev(variable),
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
    shannon = function(df_shannon, col_name, Brewer.palID="Set2") {
      if (!is(df_shannon, "data.table")) {
        stop("shannon_df needs to be a data.table")
      } else {
        # Pivot into long table
        shannon_long <- data.table::melt(data = df_shannon,
                                         measure.vars = colnames(df_shannon)[grepl("depth-", colnames(df_shannon))], 
                                         variable.name = "iters",
                                         variable.factor = FALSE,
                                         value.name = "alpha_div")
        # Corrects colnames
        colnames(shannon_long) <- c("SAMPLE-ID", "iters", "alpha_div")
        # Adds new column
        shannon_final <- base::merge(shannon_long, 
                                     self$metaData[, .SD, .SDcols = c("SAMPLE-ID", col_name)], 
                                     by = "SAMPLE-ID", 
                                     all.x = TRUE)
        
        # NAs may appear since metadata and shannon file do not have the same IDs after subsetting
        shannon_final <- na.omit(shannon_final)
        
        # Creating color palette
        unique_groups <- unique(self$metaData[[col_name]])
        chosen_palette <- RColorBrewer::brewer.pal(length(unique_groups), Brewer.palID)
        colors <- stats::setNames(chosen_palette, unique_groups)
        
        # Creating shannon plot
        return(
          shannon_final %>%
            ggplot(mapping = aes(x = base::get(col_name, shannon_final),
                                 y = alpha_div)) +
            geom_violin(width = 1.4, aes(fill = base::get(col_name, shannon_final))) +
            geom_boxplot(width = 0.1) +
            theme_bw() +
            theme(legend.position = "none",
                  text = element_text(size = 12, color = "black")) + 
            scale_fill_manual(name = "", 
                              values = colors) +
            labs(title = NULL,
                 subtitle = paste0("selected column: ", col_name),
                 x = "sample groups",
                 y = "Shannon Index")
        )
      }
    },
    composition = function(feature_rank, feature_filter = NA, col_name = NA, feature_top = 10, Brewer.palID = "RdYlBu") {
      # Agglomerate by feature_rank
      self$feature_glom(feature_rank = feature_rank, feature_filter = feature_filter)
      
      # Normalizes sample counts
      self$transform(function(x) x / sum(x))
      
      # Fetch unfiltered and filtered features
      dt <- self$countData[, (feature_rank) := self$featureData[[feature_rank]]]
      
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
      
      # returns results as list
      return(
        list(
          data = composition_final,
          palette = taxa_colors_ordered
        )
      )
    },
    ordination = function(metric, method, group_by, distmat = NULL, weighted = FALSE, normalize = TRUE, parallel = FALSE, 
                          pca.pairwise = FALSE, pca.max.explained = 80, pca.dim = c(1,2), outdir=".", cpus = 8) {
      
      if (parallel == TRUE) {
        # Uses available CPUs for %dopar%
        RcppParallel::setThreadOptions(numThreads = cpus)
      }
      
      # Normalizes counts
      if (normalize == TRUE) {
        self$transform(function(x) x / sum(x))
      }
      
      # computes distance matrix without sample rarefying
      if (is.null(distmat)) {
        # Requires rownames to contain same labels as tree
        counts <- self$countData
        counts <- as.matrix(counts)
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
        "pcoa" = pairwise.adonis(distmat, factors = self$metaData[[ group_by ]]),
        "nmds" = pairwise.anosim(distmat, grouping = self$metaData[[ group_by ]])
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
                                           Label = "p.adjusted",
                                           Y_title = "Pseudo F test statistic",
                                           plot.title = "PERMANOVA")
        # Loading score plot
        plot_list$scores_plot <- ordination_plot(df_pcs_points, 
                                                 pcs, 
                                                 pair=c("PC1", "PC2"), 
                                                 metric, 
                                                 group_by)
        
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
                                                 metric, 
                                                 group_by)
      }
      return(plot_list)
    },
    differential_feature_expression = function() {
      # place holder
    },
    correlation = function() {
      # place holder
    }
  )
)