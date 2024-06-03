ps_composition <- function(ps, tax_level, metadata.columns, Brewer.palID = "RdYlBu",taxa_n=15, excel_path=FALSE) {
  
  # Fetch OTU table and Metadata
  otu_tab <- get_otu(ps, tax_target = tax_level)
  meta_tab <- get_meta(ps)
  
  df <- as.data.frame(t(otu_tab))
  df$Taxa <- rownames(df)
  
  # Select top taxa
  df_clean <- df %>% 
    filter(!grepl("uncultured", Taxa), !grepl("metagenome", Taxa)) %>% 
    select(-Taxa)
  df_topTaxa <- df_clean %>% 
    top_taxa(n = taxa_n)
  
  # Sorts taxa by alphabetic order
  df_topTaxa$Taxa <- rownames(df_topTaxa)
  df_topTaxa.sorted <- arrange(df_topTaxa, Taxa)
  
  # Select other taxa
  row_other <- df[!(rownames(df) %in% rownames(df_topTaxa)), ] %>% select(-Taxa) %>% colSums()
  df_final <- rbind(df_topTaxa.sorted, row_other)
  rownames(df_final)[nrow(df_final)] <- "Other"
  df_final$Taxa <- rownames(df_final)
  
  # Creates palette
  taxa_colors_ordered <- stats::setNames(c(RColorBrewer::brewer.pal(length(df_final$Taxa)-1, Brewer.palID), "lightgrey"), df_final$Taxa)
  
  # Save as excel file
  if (excel_path != FALSE) {
    openxlsx::write.xlsx(df_final, paste0(excel_path, "/", tax_level, "_rel_abun.xlsx"),
                         asTable = TRUE,
                         rowNames = TRUE)
  }
  
  # Pivot longer table by specified id.vars
  df_melt <- reshape2::melt(df_final, id.vars = c("Taxa"))
  
  # Includes metadata content by user-specification
  df_melt <- add_metadata(df_melt, meta_tab, meta_col.id = "SAMPLE.ID", meta_col.add = metadata.columns)

  # Change names: HARDCODED
  #----------------------------------------------------------------------------#
  df_melt$variable <- gsub("_\\d+", "", df_melt$variable)
  #----------------------------------------------------------------------------#
  # Factors the melted dataframe by the original order of Taxa
  # Important for scale_fill_manual taxa order
  df_melt$Taxa <- factor(df_melt$Taxa, levels = df_final$Taxa)
  # composition relative abundance
  result <- list(
    df = df_melt,
    palette = taxa_colors_ordered,
    tax_rank = tax_level
  )
  return(result)
}
