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
      id_list <- data.table::copy(self$featureData[, ID])
      counts <- data.table::copy(self$countData[, ID := id_list])
      features <- data.table::copy(self$featureData)
      
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
    alpha_diversity = function(custom_div = NA, col_name, method = c("shannon", "invsimpson", "simpson"), Brewer.palID="Set2", evenness = FALSE) {
      # TO DO:
        # - Add hillR::hill_func_parti_pairwise(comm = counts, traits = metadata, q = 2)
        # - Find way to summarize hill results
      
      # OUTPUT: Plot list
      plot_list <- list(
        diversity = NULL,
        fisher_alpha = NULL
      )
      
      # Save tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData
      )
      
      # Compute diversity and other metrics if custom_div is empty
      if (is.na(custom_div)) {
        # Get matrix
        mat <- as.matrix(self$countData)
        
        # Alpha diversity based on 'method'
        div <- data.table::data.table(vegan::diversity(t(mat), index=method))
        div[, (paste(col_name)) := self$metaData[, .SD, .SDcols = c(col_name)]]
        # Adjusts for evenness
        if (evenness) div$V1 <- div$V1 / log(vegan::specnumber(div$V1)) 
        
        # Fisher alpha based on 'method'
        fish <- data.table::data.table(vegan::fisher.alpha(t(mat), index=method))
        fish[, (paste(col_name)) := self$metaData[, .SD, .SDcols = c(col_name)]]
        # Adjusts for evenness
        if (evenness) fish$V1 <- fish$V1 / log(vegan::specnumber(fish$V1)) 
        
        # get colors
        colors <- fetch_palette(self$metaData, col_name, Brewer.palID)
        
        # Create and saves plots
        plot_list$diversity <- diversity_plot(dt = div,
                                              values = "V1",
                                              col_name = col_name,
                                              palette = colors,
                                              method = method)
        plot_list$fisher_alpha <- diversity_plot(dt = fish,
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
    composition = function(feature_rank, feature_filter = NA, col_name = NA, feature_top = 10, Brewer.palID = "RdYlBu") {
      # Copies object to prevent modification of tools class components
      private$tmp_link(
        .countData = self$countData,
        .featureData = self$featureData,
        .metaData = self$metaData
      )

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
    ordination = function(metric, method, group_by, distmat = NULL, weighted = FALSE, normalize = TRUE, parallel = FALSE, 
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
                                     tree = tree)
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
                                                 metric, 
                                                 group_by)
      }
      
      # Restores tools class components
      private$tmp_restore()
      
      return(plot_list)
    },
    differential_feature_expression = function(feature_rank, sample.id, paired, paired.id, 
                                               condition.group, condition_A, condition_B, 
                                               feature_filter = NA, feature_top = 20, normalize = TRUE, cpus = 8) {
      # Final output
      plot_list <- list(
        data = NULL,
        boxplot = NULL,
        barplot = NULL,
        tile_plot = NULL,
        rel_abun = NULL,
        volcano_plot = NULL,
        pvalues = NULL,
        volcano = NULL
      )
      
      # Copies object to prevent modification of tools class components
      counts <- data.table::copy(self$countData)
      metadata <- data.table::copy(self$metaData)
      features <- data.table::copy(self$featureData)
      
      #------#
      # Main #
      #------#
      
      # Subset samples by conditions
      # self$sample_subset(self$metaData[[condition.group]] %in% c(condition_A, condition_B))
      
      if (normalize) {
        self$transform(function(x) x / sum(x))
      }
      
      # Agglomerate taxa by feature rank and filter unwanted taxa
      self$feature_glom(feature_rank = feature_rank, feature_filter = feature_filter)
      
      # Creates long table of relative abundance
      dt <- self$countData[, (feature_rank) := self$featureData[[feature_rank]]]
      stats_dt <- base::merge(data.table::melt(dt,
                                               measure.vars = colnames(dt)[!grepl(feature_rank, colnames(dt))],
                                               variable.name = sample.id, 
                                               value.name = "values"),
                              self$metaData[, .SD, .SDcols = c(sample.id, condition.group)],
                              by = sample.id)
      
      # Create row_sums
      dt[, row_sum := rowSums(.SD), .SDcols = !c(feature_rank)]
      
      # Check how many features to select (depended if volcano is desired)
      if (!is.na(feature_top)) {
        feature_top <- feature_top
      } else {
        feature_top <- nrow(self$featureData)
      }
      
      # Orders by row_sum in descending order
      self$countData <- data.table::setorder(dt, -row_sum)[1:feature_top, .SD, .SDcols = !c("row_sum")]
      features <- self$countData[[ feature_rank ]]
      self$countData <- self$countData[, .SD, .SDcols = !c(feature_rank)]
      
      # Log2 transform taxa
      self$transform(log2)
      
      # Subset by top features
      stats_dt <- stats_dt[stats_dt[[feature_rank]] %in% features]
      dt <- self$countData[, (feature_rank) := features]
      
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
        # Save data
        plot_list$data <- DFE$data
        
      } else {
        stop("paired can only be TRUE or FALSE, check your input.")
      }
      
      # First merge 
      add_columns <- unique(taxa$metaData[, .(sample.id, paired.id)])
      
      merged_data <- base::merge(
        stats_dt,
        add_columns,
        by = sample.id,
        all.x = TRUE
      )
      
      # Subset merged data
      subset_merged <- merged_data[, .(patient.id, feature_rank, values)]
      colnames(subset_merged) <- c("SAMPLE-ID", "Genus", "values")
      
      # Second merge
      final_merge <- base::merge(
        x = DFE,
        y = subset_merged,
        by = c("SAMPLE-ID", "Genus"),
        all.x = TRUE
      )
      
      grouped_dt <- final_merge %>% 
        dplyr::group_by(`SAMPLE-ID`, Genus) %>% 
        dplyr::summarise(mean_values = mean(values, na.rum = TRUE),
                         diff_1 = mean(diff_1, na.rm = TRUE)) %>% 
        dplyr::ungroup()
      
      test_heatmap <- grouped_dt %>% 
        ggplot(mapping = aes(x = `SAMPLE-ID`,
                             y = Genus)) +
        geom_point(aes(size = mean_values, fill = diff_1), shape = 21) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
              axis.text.y = element_text(size=12),
              axis.text = element_text(size=12),
              text = element_text(size=12),
              legend.text = element_text(size=12),
              legend.title = element_text(size=14),
              strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
        scale_y_discrete(limits = rev(levels(as.factor(grouped_dt[["Genus"]])))) +
        scale_fill_gradient2(name = paste0("log2( C5 / C1 )"),
                             low = "blue",
                             mid = "white",
                             high = "red",
                             na.value = "grey80") +
        scale_size_continuous(name = "Mean Rel. Abun. (%)", labels = scales::label_number(accuracy = 0.01)) +
        labs(x = NULL, 
             y = NULL)
      
      # Generate heatmap plot with df_diff data
      if (paired == TRUE) {
        # Generating heatmap plot based on paired boolean
        n_diff_columns <-  sum(grepl("^diff_", colnames(DFE)))
        
        # Generate heatmap plot with df_diff data
        heatmap_plot <- DFE %>% 
          ggplot(mapping = aes(x = base::get(sample.id, DFE),
                               y = base::get(feature_rank, DFE)))
        
        # If there is only one column uses default settings
        if (n_diff_columns == 1) {
          heatmap_plot <- heatmap_plot +
            geom_tile(aes(fill = diff_1))
        } else {
          # Adds geom_tile for number of diff_columns
          for (i in 1:n_diff_columns) {
            if (i == 1) {
              heatmap_plot <- heatmap_plot +
                geom_tile(aes(fill = !!sym(paste0("diff_", i))), width = 0.45)
            } else {
              heatmap_plot <- heatmap_plot +
                geom_tile(aes(fill = !!sym(paste0("diff_", i))), width = 0.45,
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
                axis.title.y = element_blank(),
                strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
          scale_fill_gradient2(name = paste0("log2( A / B )"),
                               low = "blue",
                               mid = "white",
                               high = "red",
                               na.value = "grey80") +
          scale_y_discrete(limits = rev(levels(as.factor(DFE[[ feature_rank ]])))) +
          labs(x = NULL, 
               y = "Features")  
        
        
        # Creates boxplot from relative abundances
        plot_list$rel_abun <- stats_dt %>% 
          ggplot(mapping = aes(x = values,
                               y = .data[[ feature_rank ]])) +
          facet_wrap(~.data[[condition.group]], ncol = length(condition_A) + length(condition_B)) +
          geom_boxplot() +
          theme_bw() +
          theme(text=element_text(size=12),
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                panel.spacing.x = unit(1, "lines")) +
          scale_x_continuous(trans = scales::log_trans()) +
          labs(x = "Log10( Rel. Abun. )")
      } else {
        # Store pvalues and volcano dt
        plot_list$volcano <- DFE$volcano
        plot_list$pvalues <- DFE$pvalues
        # Generating heatmap plot based on paired boolean
        n_diff_columns <-  sum(grepl("^diff_", colnames(DFE$data)))
        # Creates boxplot and barplot for unpaired samples
        for (k in c("boxplot", "barplot")) {
          plot_list[[k]] <- patchwork::wrap_plots(
            lapply(1:n_diff_columns,
                   function(i) fold_plot(dt = DFE$data, 
                                         X = paste0("diff_", i), 
                                         Y = feature_rank,
                                         title = paste0("Log2 ( ", condition_A[i], " / ", condition_B[i], " )"), 
                                         method = k, 
                                         taxa_labels = i == 1, 
                                         pvalues = DFE$pvalues)),
            ncol = n_diff_columns,
            nrow = 1)
        }
      }
      
      # Restores tools class components
      self$countData <- counts
      self$featureData <- features
      self$metaData <- metadata
      
      return(plot_list)
    },
    triplot = function() {
      # place holder for general triplot function, takes any type of model
    },
    correlation = function() {
      # place holder
    },
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