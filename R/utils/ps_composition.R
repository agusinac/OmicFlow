ps_composition <- function(ps, tax_level, metadata.columns, Brewer.palID = "RdYlBu",taxa_n=15, excel_path=FALSE) {
  
  # Fetch OTU table and Metadata
  otu_tab <- get_otu(ps, tax_target = tax_level, top_n = FALSE, filter_taxa = FALSE)
  meta_tab <- get_meta(ps)
  
  df <- as.data.frame(t(otu_tab))
  df$Taxa <- rownames(df)

  # Select top taxa
  df_clean <- df %>% 
    filter(!grepl("uncultured|metagenome", Taxa)) %>% 
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
  df_taxa_len <- length(df_final$Taxa)
  if (Brewer.palID == FALSE) {
    chosen_palette <- viridis::viridis(df_taxa_len - 1)
  } else if (df_taxa_len-1 <= 15 & df_taxa_len-1 > 10) {
    chosen_palette <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                        "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                        "#920000","#924900","#db6d00","#24ff24","#ffff6d")[1:df_taxa_len-1]
  } else {
    chosen_palette <- RColorBrewer::brewer.pal(df_taxa_len-1, Brewer.palID)
  }
  taxa_colors_ordered <- stats::setNames(c(chosen_palette, "lightgrey"), df_final$Taxa)
  
  # Save as excel file
  if (excel_path != FALSE) {
    openxlsx::write.xlsx(df_final, paste0(excel_path, "/", tax_level, "_rel_abun.xlsx"),
                         asTable = TRUE,
                         rowNames = TRUE)
  }
  
  # Pivot longer table by specified id.vars
  df_melt <- reshape2::melt(df_final, id.vars = c("Taxa"))
  
  # Includes metadata content by user-specification
  for (variable in metadata.columns) {
    df_melt <- df_melt %>% 
      rowwise() %>% 
      mutate(
        !!variable := meta_tab[[ {{ variable }} ]][stringr::str_detect(as.character(meta_tab[["SAMPLE.ID"]]), as.character(variable))]
      )
  }

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
