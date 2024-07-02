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
#
#
ps_fold_plot <- function(ps, taxa_n = 20, taxa_rank = "Genus", col_id = "PATIENT.ID", col_group = "RANKSTAT_treatment", condition_A, condition_B, method = "paired") {
  
  # Final output
  plot_list <- list(
    data = NULL,
    boxplot = NULL,
    barplot = NULL,
    tile_plot = NULL,
    rel_abun = NULL
  )
  
  paired_fold <- function(df_final, condition_A, condition_B) {
    for (i in 1:length(condition_A)) {
      name <- paste0("diff_",i)
      stats_name <- paste0("p.value_diff_",i)
      df_diff <- df_final %>% 
        group_by(sample.id, Taxa) %>% 
        summarise(!!name := value[group == condition_A[i]] - value[group == condition_B[i]])
      
      df_final <- df_final %>% 
        left_join(df_diff, by = c("Taxa", "sample.id"))
    }
    
    result <- list(
      df = df_final,
      diff_n = sum(grepl("^diff_", names(df_final)))
    )
    
    return(result)
  }
  
  unpaired_fold <- function(df_final, condition_A, condition_B) {
    # Creates a memory efficient data frame
    data.table::setDT(df_final)
    dt_A <- df_final[group %in% condition_A, .(Taxa, value, group)]
    dt_B <- df_final[group %in% condition_B, .(Taxa, value, group)]
    
    # Initialize an empty list to store the results
    result_list <- vector("list", length(condition_A))
    for (i in 1:length(condition_A)) {
      name <- paste0("diff_", i)
      
      # Filter dt_A and dt_B for the current conditions
      dt_A_i <- dt_A[group == condition_A[i]]
      dt_B_i <- dt_B[group == condition_B[i]]
      
      # Perform the join and subtraction
      dt_diff <- dt_A_i[dt_B_i, on = .(Taxa), allow.cartesian=TRUE]
      dt_diff[, (name) := value - i.value]
      dt_diff[, (paste0("p.value_diff_", i)) := wilcox.test(value - i.value, correct = TRUE)$p.value]
      
      # Store the result in the list
      result_list[[i]] <- dt_diff
    }
    # Combine the results into a single data.table
    df_joined <- data.table::rbindlist(result_list, fill = TRUE)
    
    result <- list(
      df = df_joined,
      diff_n = sum(grepl("^diff_", names(df_joined)))
    )
    
    return(result)
  }
  
  fold_plot <- function(df, X, title, method, taxa_labels = FALSE) {
    if (method == "barplot") {
      plt <- df %>% 
        aggregate(formula(paste(X, " ~ Taxa")), sum) %>% 
        ggplot(mapping = aes(x = base::get(X),
                             y = Taxa,
                             fill = base::get(X))) +
        geom_bar(stat = "identity")
    } else if (method == "boxplot") {
      plt <- df %>% 
        ggplot(mapping = aes(x = base::get(X),
                             y = Taxa)) +
        geom_boxplot()
    }
    
    plt <- plt + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
            axis.text.y = element_text(size=12),
            axis.text = element_text(size=12),
            text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            legend.position = "none",
            axis.title.y = element_blank(),
            strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))
    if (!taxa_labels) {
      plt <- plt + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    }
    plt <- plt +
      scale_fill_gradient2(name = paste0("log2( A / B )"),
                           low = "blue",
                           mid = "white",
                           high = "red",
                           na.value = "grey80") +
      scale_y_discrete(limits = rev(levels(as.factor(df$Taxa)))) +
      labs(x = NULL,
           y = "Taxa") +
      ggtitle(title)
    
    return(plt)
  }
  #############
  # MAIN CODE #
  #############
  
  # Fetch otus and filter out less representative Genera
  otu_tab <- get_otu(ps, tax_target = taxa_rank, top_n = taxa_n, filter_taxa = TRUE) %>% t()
  meta_tab <- get_meta(ps)
  
  # log2 scale
  otu_tab.trans <- otu_tab %>%
    as.matrix() %>% 
    as("sparseMatrix") %>% 
    log2() %>% 
    as.matrix()
  
  # replace infinite by zero
  otu_tab.trans[!is.finite(otu_tab.trans)] <- 0
  
  # Creates dataframe & orders taxa abundances
  df <- as.data.frame(otu_tab.trans)
  df$Taxa <- rownames(df) 
  df$Total <- rowSums(df[, -ncol(df)])
  df <- df[order(-df$Total, decreasing=TRUE), ]
  
  # Pivot longer table by specified id.vars
  df_melt <- reshape2::melt(df, id.vars = c("Taxa", "Total"))
  
  
  # Creates long table with or without paired data
  df_final <- df_melt %>% 
    rowwise() %>% 
    mutate(
      # Collects sample group names from metadata 
      group = meta_tab[[ {{ col_group }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))],
      # Metadata should have a column of identical sample
      sample.id = meta_tab[[ {{ col_id }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))]
    )
  
  if (is.null(df_final$group)) {stop(paste("Col_group isnt found in metadata, check spelling!"))}
  
  # Removes unpaired samples
  if (method == "paired") {
    if (is.null(df_final$sample.id)) {stop(paste("col.id is not found in metadata, check spelling!"))}
    df_final <- df_final %>% 
      group_by(sample.id) %>% 
      filter(any(group %in% {{condition_A}}) & any(group %in% {{condition_B}}))
  }
  
  
  # Computes log2(A) - log2(B)
  # Supports multiple inputs for A and B.
  # For example A = T1, T2 and B = H1, H2
  
  # Efficient memory data.table
  method_data <- switch(
    method,
    "paired" = paired_fold(df_final, condition_A, condition_B),
    "unpaired" = unpaired_fold(df_final, condition_A, condition_B)
  )
  
  # To be used for paired or unpaired
  df_final <- method_data$df %>% as.data.frame()
  n_diff_columns <- method_data$diff_n
  
  # Save data
  plot_list$data <- df_final
  
  # Generate heatmap plot with df_diff data
  if (method == "paired") {
    # Generate heatmap plot with df_diff data
    heatmap_plot <- df_final %>% 
      ggplot(mapping = aes(x = as.factor(sample.id),
                           y = Taxa))
    
    # If there is only one column uses default settings
    if (n_diff_columns == 1) {
      heatmap_plot <- heatmap_plot +
        geom_tile(aes(fill = !!sym(paste0("diff_", i))))
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
      scale_y_discrete(limits = rev(levels(as.factor(df_final$Taxa)))) +
      labs(x = NULL, 
           y = "Taxa")
    
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
    
    # Creates boxplot from relative abundances
    plot_list$rel_abun <- stats_final %>% 
      ggplot(mapping = aes(x = value,
                           y = Taxa)) +
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
                                     taxa_labels = i == 1)),
        ncol = n_diff_columns,
        nrow = 1)
    }
    
  }
  return(plot_list) 
}