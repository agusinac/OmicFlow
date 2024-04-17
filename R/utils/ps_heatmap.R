ps_fold_heatmap <- function(ps, threshold = 0.2, taxa_rank = "Genus", negative_ct="healthy", positive_ct="tumor") {
  # Fetch otus and filter out less representative Genera
  otu_tab <- ps2otu(ps, tax_target = taxa_rank) %>% t()
  otu_tab.filt <- otu_tab[rowSums(otu_tab) > threshold, ]
  
  # log2 scale
  otu_tab.trans <- otu_tab.filt %>%
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
  
  # Removes uncultured genera, cleans sample names
  df_final <- df_melt %>% 
    filter(!grepl("uncultured", Genera)) %>% 
    mutate(
      # Specific for data, atm hardcoded but should be automated!
      group = ifelse(base::grepl("T_", variable), "tumor",
                     ifelse(base::grepl("H_", variable), "healthy", NA)),
      sample.id = sub("(-?\\w{1}_\\d+)$", "" ,variable)
    ) %>% 
    group_by(sample.id) %>% 
    filter(any(group == negative_ct) & any(group == positive_ct))
  
  # Computes log2(tumor) - log2(healthy)
  df_diff <- df_final %>% 
    group_by(sample.id, Genera) %>% 
    summarise(diff = value[group == positive_ct] - value[group == negative_ct])
  
  # Generate heatmap plot with df_diff data
  heatmap_plot <- df_final %>% 
    left_join(df_diff, by = c("Taxa", "sample.id")) %>% 
    ggplot(mapping = aes(x = as.factor(sample.id),
                         y = Genera,
                         fill = diff)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
    scale_fill_gradient2(name = "log2(tumor / healthy) fold",
                         low = "blue",
                         mid = "white",
                         high = "red") +
    scale_y_discrete(limits = rev(levels(as.factor(df_final$Genera)))) +
    labs(x = "Samples", 
         y = "Taxa")
  
  # Fetch otu table for boxplot and reshapes into long table
  otu_tab.tumor <- as.data.frame(otu_tab.filt)
  otu_tab.tumor$Genera <- rownames(otu_tab.filt)
  otu_tab.tumor_ <- reshape2::melt(otu_tab.tumor, id.vars = c("Taxa"))
  
  rel_abun_plot <- otu_tab.tumor_ %>% 
    ggplot(mapping = aes(x = value,
                         y = Genera)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_continuous(breaks=seq(0, 1, 0.5)) +
    scale_y_discrete(limits = rev(levels(as.factor(df_final$Genera)))) +
    labs(x = "Relative Abundance")
  
  comb_plot <- (heatmap_plot + rel_abun_plot) +
    plot_layout(guides = "collect",
                axes = "collect",
                widths = c(5,2)) +
    plot_annotation(title = "Microbiome of Breast Carcinoma patients",
                    subtitle = "")
  
  return(comb_plot)
}