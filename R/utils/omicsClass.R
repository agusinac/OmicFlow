### TO DO:
# - Add functions:  -> composition
# - Add functions:  -> ordinations
# - Add functions:  -> fold plot
# - Add functions:  -> Linear regressin (incl. RDA, CCA etc.)

# Required libraries to be loaded:
library("foreach")
library("tidyr")
library("ggplot2")
library("magrittr")

# R files that are still required
source("automated-omics-analysis/R/utils/pairwise.adonis.R")
source("automated-omics-analysis/R/utils/pairwise.anosim.R")
source("automated-omics-analysis/R/utils/parse_commandline.R")
source("automated-omics-analysis/R/utils/composition_plot.R")

# Additional functions removed from nested methods
source("automated-omics-analysis/R/utils/ordination_plot.R")
source("automated-omics-analysis/R/utils/stats_plot.R")
source("automated-omics-analysis/R/utils/paired_fold.R")
source("automated-omics-analysis/R/utils/unpaired_fold.R")
source("automated-omics-analysis/R/utils/DFE_plot.R")


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
    },
    feature_subset = function(...) {
      rows_to_keep <- self$featureData[, ...]
      self$featureData <- self$featureData[rows_to_keep, ]
      self$countData <- self$countData[rows_to_keep, ]
      self$removeZeros()
      invisible(self)
    },
    sample_subset = function(...) {
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
    composition = function(feature_rank, feature_filter = NA, col_name, feature_top = 10, Brewer.palID = "RdYlBu"
) {
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
      final_dt <- rbind(clean_dt[1:feature_top][, .SD, .SDcols = !c("row_sum")], 
                        clean_dt[(feature_top+1):nrow(clean_dt)][, lapply(.SD, function(x) sum(x)), 
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
      composition_final <- base::merge(final_long,
                                       self$metaData[, .SD, .SDcols = c("SAMPLE-ID", col_name)],
                                       by = "SAMPLE-ID",
                                       all.x = TRUE)
      
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
                          pca.pairwise = FALSE, pca.max.explained = 80, pca.dim = c(1,2), outdir=".") {

      if (parallel == TRUE) {
        # Uses available CPUs for %dopar%
        doParallel::registerDoParallel(cores = parallel::detectCores())
      }
      
      # Normalizes counts
      if (normalize == TRUE) {
        self$transform(function(x) x / sum(x))
      }
      
      # computes distance matrix without sample rarefying
      if (is.null(distmat)) {
        # Creates empty matrix
        n <- ncol(self$countData)
        sample_names <- colnames(self$countData)
        dist <- matrix(0, nrow = n, ncol = n, 
                       dimnames = list(sample_names, sample_names))
        
        pairs <- utils::combn(x = sample_names, 
                              m = 2, 
                              simplify = FALSE)
        
        # Duplicates matrix to prevent over-use with $ and parellel compatibility.
        count_mat <- self$countData
        tree <- self$treeData
        
        # Performs sample by sample computation. Uses CPUs if enabled by parallel
        distlist <- foreach::foreach(pair = pairs) %dopar% {
          a <- count_mat[[ pair[1] ]]
          b <- count_mat[[ pair[2] ]]
          
          # weighted / unweighted
          if (weighted == FALSE) {
            a[a > 0] <- 1
            b[b > 0] <- 1
          }
          
          # Switch statement for different (dis)similarity metrics
          d <- switch(
            metric,
            "unifrac" = if (is(tree, "phylo")) {
              abdiv::weighted_unifrac(x = a,
                                      y = b,
                                      tree = tree)
            } else {
              stop("tree data is missing")
            },
            "bray" = abdiv::bray_curtis(x = a, y = b),
            "jaccard" = abdiv::jaccard(x = a, y = b)
          )
        }
        
        # convert parts to matrix
        dist[lower.tri(dist)] <- as.numeric(distlist)
        distmat <- t(dist)
        distmat[lower.tri(distmat)] <- as.numeric(distlist)
      }
      distmat <- as.dist(distmat)
      
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


# Testing class:
setwd(paste0(getwd(), "/automated-omics-analysis/data/"))
test <- metataxonomics$new(metaData = "metadata.tsv",
                           biomData = "biom_with_taxonomy.biom",
                           treeData = "rooted_tree.newick")

test$feature_subset(Domain == "Bacteria")
# differential feature expression 

differential_feature_expression = function() {
  # Required parameters
  feature_rank <- "Genus"
  sample.id <- "SAMPLE-ID"
  paired <- TRUE
  paired.id <- "PATIENT-ID"
  feature_filter <- NA
  feature_top <- 20
  condition.group <- "RANKSTAT_treatment"
  condition_A <- c("T")
  condition_B <- c("H")
  normalize <- TRUE
  
  # Final output
  plot_list <- list(
    data = NULL,
    boxplot = NULL,
    barplot = NULL,
    tile_plot = NULL,
    rel_abun = NULL,
    volcano_plot = NULL
  )
  
  #------#
  # Main #
  #------#
  
  if (normalize) {
    test$transform(function(x) x / sum(x))
  }
  
  # Agglomerate taxa by feature rank and filter unwanted taxa
  test$feature_glom(feature_rank = feature_rank, feature_filter = feature_filter)
  
  # keeps copy of relative abundance
  counts <- test$countData
  stats_rel <- counts[, (feature_rank) := test$featureData[[feature_rank]]]
  
  # Log2 transform taxa
  test$transform(log2)
  
  # Sort and select top features
  dt <- test$countData[, (feature_rank) := test$featureData[[feature_rank]]]
  
  # Create row_sums
  dt[, row_sum := rowSums(.SD), .SDcols = !c(feature_rank)]
  
  # Orders by row_sum in descending order
  dt <- data.table::setorder(dt, -row_sum)[1:feature_top, .SD, .SDcols = !c("row_sum")]
  
  # Copying metadata to prevent in place modifications
  metadata <- test$metaData

  # Compute 2-fold expression based on (un)paired samples
  # Computes on equation oflog2(A) - log2(B)
  # Supports multiple inputs for A and B.
  # For example A = T1, T2 and B = H1, H2
  if (paired == TRUE) {
    # sorting of metadata
    condition.labels <- data.table::setorderv(metadata, 
                                              cols = c(sample.id, paired.id, condition.group))[[ condition.group ]]
    # paired samples
    DFE <- paired_fold(dt = dt,
                       sample.id = sample.id,
                       paired.id = paired.id,
                       condition_A = condition_A,
                       condition_B = condition_B,
                       unique.id = unique(test$metaData[[ paired.id ]]),
                       condition_labels = condition.labels,
                       feature_rank = feature_rank)
  } else if (paired == FALSE) {
    # sorting of metadata
    condition.labels <- data.table::setorderv(metadata,
                                              cols = c(sample.id, condition.group))[[ condition.group ]]
    # unpaired samples
    DFE <- unpaired_fold(dt = dt,
                         sample.id = sample.id,
                         condition_A = condition_A,
                         condition_B = condition_B,
                         condition_labels = condition.labels,
                         feature_rank = feature_rank)
  } else {
    stop("paired can only be TRUE or FALSE, check your input.")
  }

  # Generating heatmap plot based on paired boolean
  n_diff_columns <-  sum(grepl("^diff_", colnames(DFE)))
  
  # Generate heatmap plot with df_diff data
  if (paired == TRUE) {
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
    stats_dt <- data.table::melt(stats_rel,
                                 measure.vars = colnames(stats_rel)[!grepl(feature_rank, colnames(stats_rel))],
                                 variable.name = "samples", 
                                 value.name = "values")
    
    
    plot_list$rel_abun <- stats_dt %>% 
      ggplot(mapping = aes(x = base::get(value, stats_dt),
                           y = base::get(feature_rank, stats_dt))) +
      geom_boxplot() +
      facet_wrap(~group, ncol = length(condition_A) + length(condition_B)) +
      theme_bw() +
      theme(text=element_text(size=12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.spacing.x = unit(1, "lines")) +
      scale_x_continuous(trans = scales::log_trans()) +
      scale_y_discrete(limits = rev(levels(as.factor(df_final$Taxa)))) +
      labs(x = "Log10( Rel. Abun. )")
  
  
    
    
  
}


# Creates Differential 2-fold change expression from phyloseq object
# 
# Works on both paired and non-paired data
# 
# Returns a list of following items:
#   - data frame
#   - Boxplot of fold expressions for different groups
#   - Barplot of Summed fold expressions 
#   - Heatmap plot for paired samples
#   - Relative abundance distribution of paired samples

ps_fold_plot <- function(ps, taxa_n = 20, taxa_rank = "Genus", col_id = "PATIENT.ID", col_group = "RANKSTAT_treatment", condition_A, condition_B, method = "paired", stat_test = FALSE) {
  
  #############
  # MAIN CODE #
  #############
# To be converted into data.table operations
  if (stat_test) {
    # Compute wilcox significance between taxa for each condition_A
    taxa_groups <- unique(df_final$Taxa)
    taxa_pvalues <- matrix(list(), ncol = length(condition_A), nrow = length(taxa_groups))
    
    for (i in 1:length(condition_A)) {
      for (j in 1:length(taxa_groups)) {
        sub.df <- filter(df_final, Taxa == taxa_groups[j] & group == condition_A[i])
        taxa_pvalues[j, i] <- stats::wilcox.test(sub.df$value, sub.df$i.value, correct = TRUE)$p.value
      }
    }
    # Data wrangling and merging to df_final
    rownames(taxa_pvalues) <- taxa_groups
    
    if (ncol(taxa_pvalues) == 1) {
      taxa_pvalues.df <- t(as.data.frame(taxa_pvalues[sort(rownames(taxa_pvalues)), ]))
    } else {
      taxa_pvalues.df <- as.data.frame(taxa_pvalues[sort(rownames(taxa_pvalues)), ])
    }
    colnames(taxa_pvalues.df) <- condition_A
    df_final <- base::merge(df_final, taxa_pvalues,
                            by.x = "Taxa", by.y = "row.names", all.x = TRUE)
    df_final[, c("V1", "V2")]
  }

  
  
  
  
    
    # Fetch otu table for boxplot and reshapes into long table
    stats_tab <- as.data.frame(otu_tab)
    stats_tab$Taxa = rownames(stats_tab)
    
    # Pivot longer
    stats_melt <- reshape2::melt(stats_tab, id.vars = c("Taxa"))
    
    # Validate numeric zero's instead of NAs
    stats_melt$value <- as.numeric(stats_melt$value)
    stats_melt$value[is.na(stats_melt$value)] <- 0
    
    stats_final <- stats_melt %>% 
      rowwise() %>% 
      mutate(
        # Collects sample group names from metadata 
        sample.id = meta_tab[[ {{ col_id }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))],
        group = meta_tab[[ {{ col_group }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))]
      ) %>% 
      group_by(sample.id) %>% 
      filter(any(group %in% {{condition_A}}) & any(group %in% {{condition_B}}))
    

    # Uses paired samples also to view box/barplot of whole groups
    for (k in c("boxplot", "barplot")) {
      plot_list[[k]] <- patchwork::wrap_plots(
        lapply(1:n_diff_columns,
               function(i) fold_plot(df = df_final, 
                                     X = paste0("diff_", i), 
                                     title = paste0("Log2 ( ", condition_A[i], " / ", condition_B[i], " )"), 
                                     method = k, 
                                     taxa_labels = i == 1)),
        ncol = n_diff_columns,
        nrow = 1)
    }
    
  } else if (method == "unpaired") {
    # Graphs are duplicate for each method, Save each graph into a list and then arrange graph at the end!
    for (k in c("boxplot", "barplot")) {
      plot_list[[k]] <- patchwork::wrap_plots(
        lapply(1:n_diff_columns,
               function(i) fold_plot(df = df_final, 
                                     X = paste0("diff_", i), 
                                     title = paste0("Log2 ( ", condition_A[i], " / ", condition_B[i], " )"), 
                                     method = k, 
                                     taxa_labels = i == 1, 
                                     pvalues = taxa_pvalues.df,
                                     pvalues.col = condition_A[i])),
        ncol = n_diff_columns,
        nrow = 1)
    }
  }
  return(plot_list) 
}
  

# measure memory of object:
pryr::object_size(test)
pryr::object_size(ps)



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
      # check countData !
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



# Create Additional omics classes:
# - Transcriptomics
# - Metagenomics
# - Metabolomics
# - Proteomics

