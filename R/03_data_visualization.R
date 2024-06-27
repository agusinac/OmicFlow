################
### RANKSTAT ###
################

RANKSTAT_ncol <- length(RANKSTAT_data)

### Classified ASVs per taxonomic rank
#------------------------------------------------------------------------------#
proportion_classified <- ps_rankstat(ps_abs) +
  plot_annotation(title = "Number of ASVs classified at each rank") +
  plot_layout(guides = "collect",
              axis_titles = "collect")

ggsave(
  filename = "automated-omics-analysis/results/03_rankstat.png",
  plot = proportion_classified,
  limitsize = FALSE
)
#------------------------------------------------------------------------------#
# Microbiome composition by all samples


nrow <- length(taxa_names)

if (RANKSTAT_ncol > 0) {

  composition_plots <- matrix(list(), RANKSTAT_ncol, nrow)
  
  for (i in 1:RANKSTAT_ncol) {
    col_name <- colnames(RANKSTAT_data)[i]
    
    # Microbiome composition by all samples
    for (j in 1:nrow) {
      # Creates composition long table
      res <- ps_composition(ps = ps_rel,
                            tax_level = taxa_names[j],
                            metadata.columns = col_name,
                            taxa_n = 10)
      
      # Creates composition ggplot as list
      composition_plots[[i, j]] <- composition_plot(df = res$df,
                                                    palette = res$palette,
                                                    tax_level = res$tax_rank,
                                                    title_name = col_name)
    }
  }
  #Saves plot list as png file
  ggsave(filename = "automated-omics-analysis/results/03_composition_combinations.png",
         plot = patchwork::wrap_plots(composition_plots,
                                      ncol = RANKSTAT_ncol,
                                      nrow = nrow),
         width = 30,
         height = 18,
         limitsize = FALSE,
         scaling = 1)
}

###################
### CORRELATION ###
###################

CORRELATION_ncol <- length(CORRELATION_data)
cor_method <- "spearman"

if (CORRELATION_ncol > 0) {
  correlation_heatmap_plt <- list()
  
  for (j in 1:nrow) {
    # Get ASV table and compute correlation matrix
    X <- get_otu(ps, tax_target = taxa_names[j], top_n = 30) %>% as.matrix()
    cor_mat <- cor(X, CORRELATION_data, method = cor_method)
    colnames(cor_mat) <- sub("CORRELATION_", "", colnames(CORRELATION_data)) 
    
    # Subset phyloseq by taxa in correlation matrix
    ps.glom <- ps %>% 
      phyloseq::tax_glom(taxrank = taxa_names[j]) %>%
      removeZeros()
    
    ps.filt <- phyloseq::prune_taxa(phyloseq::tax_table(ps.glom)[, taxa_names[j]] %in% rownames(cor_mat),
                                    x = ps.glom)
    
    # Construct matrix from tree
    tree_mat <- ps_Tree2Matrix(ps.filt, tax_level = taxa_names[j])
    
    # Store plot in list
    correlation_heatmap_plt[[j]] <- ComplexHeatmap::Heatmap(matrix = cor_mat,
                                                            name = "Correlation Heatmap",
                                                            heatmap_legend_param = list(
                                                              title = cor_method, 
                                                              at = c(-1,-0.5,0,0.5,1),
                                                              legend_height = unit(4, "cm"),
                                                              direction = "vertical",
                                                              title_position = "topcenter"
                                                            ),
                                                            color = colorRampPalette(c("blue", "white", "red"))(100),
                                                            clustering_distance_rows = dist(tree_mat),
                                                            rect_gp = grid::gpar(col = "black",
                                                                                 lwd = 1),
                                                            row_title = taxa_names[j],
                                                            row_names_max_width = unit(6, "cm"),
                                                            column_names_gp = grid::gpar(fontsize = 6),
                                                            row_names_gp = grid::gpar(fontsize = 6),
                                                            width = unit(0.5, "cm"),
                                                            height = unit(9, "cm")
    )
  }
}


#####################
### PAIREDGROUPBY ###
#####################

PAIREDGROUPBY_ncol <- length(PAIREDGROUPBY_data)

if (PAIREDGROUPBY_ncol > 0) {

  heatmap_plots <- matrix(list(), PAIREDGROUPBY_ncol, nrow)
  
  for (i in 1:PAIREDGROUPBY_ncol) {
    # Subsets column name and data
    col_name <- colnames(PAIREDGROUPBY_data)[i]
    col_data <- PAIREDGROUPBY_data[, i]
    
    # Saves conditions A and B as list
    A <- sort(col_data[grepl("A", col_data)], decreasing = FALSE)
    B <- sort(col_data[grepl("B", col_data)], decreasing = FALSE)
    
    for (j in 1:nrow) {
      # log2 (A / B) heatmap for different taxa levels as list
      heatmap_plots[[i, j]] <- ps_heatmap(ps = ps_rel,
                                          taxa_rank = taxa_names[j],
                                          col_id = "PATIENT.ID",
                                          col_group = col_name,
                                          condition_A = A,
                                          condition_B = B
      ) + plot_annotation(subtitle = paste0("Taxa rank = ", taxa_names[j]))
    }
  }
  # Saves p list in directory as pdf
  # `gridExtra::marrangeGrob` does not work with "patchwork"
  #Saves plot list as pdf file
  ggsave(filename = "automated-omics-analysis/results/03_2fold_heatmap.png",
         plot = patchwork::wrap_plots(heatmap_plots,
                                      ncol = PAIREDGROUPBY_ncol,
                                      nrow = nrow),
         width = 15,
         height = 20,
         limitsize = FALSE,
         scaling = 1)
}
