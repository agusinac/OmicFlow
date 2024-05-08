ps_composition <- function(ps, tax_level, metadata.columns, taxa_n=15, excel_path=FALSE) {

  # Fetch OTU table and Metadata
  otu_tab <- get_otu(ps, tax_target = tax_level)
  meta_tab <- get_meta(ps)
  
  df <- as.data.frame(t(otu_tab))
  df$Taxa <- rownames(df)
  
  # Select top taxa
  df_topTaxa <- df %>% 
    filter(!grepl("uncultured", Taxa), !grepl("metagenome", Taxa)) %>% 
    select(-Taxa) %>% 
    top_taxa(n = taxa_n)
  
  df_topTaxa$Taxa <- rownames(df_topTaxa)
  
  # Save as excel file
  if (excel_path != FALSE) {
    openxlsx::write.xlsx(df_topTaxa, paste0(excel_path, "/", tax_level, "_rel_abun.xlsx"),
                         asTable = TRUE,
                         rowNames = TRUE)
  }
  
  # Pivot longer table by specified id.vars
  df_melt <- reshape2::melt(df_topTaxa, id.vars = c("Taxa"))
  
  # Includes metadata content by user-specification
  for (variable in metadata.columns) {
    df_melt <- df_melt %>% 
      rowwise() %>% 
      mutate(
        !!variable := meta_tab[[ {{ variable }} ]][str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))]
      )
  }
  
  # composition of absolute counts per sample
  #comp_abs <- sampleSums_plot(ps)
  
  # composition relative abundance
  comp_rel <- composition_plot(df_melt, tax_level)
  return(comp_rel)
}