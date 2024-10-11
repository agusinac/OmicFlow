rda_select <- function(otu_tab.log, meta_table) {
  # Collect RANKSTAT_ column names
  rankstat_meta_tab <- meta_table[, grepl("RANKSTAT_", names(meta_table))] 
  
  # Checks if table contains at least 2 variables
  if (!is.null(dim(rankstat_meta_tab))) {
    # Creates model against 1 (simplistic model)
    mod0 <- vegan::rda(otu_tab.log ~ 1, 
                       data = rankstat_meta_tab,
                       scale = FALSE, 
                       na.action = na.fail, 
                       subset = NULL
    )
    
    # Full model against all variables
    mod1 <- vegan::rda(otu_tab.log ~ ., 
                       data = rankstat_meta_tab,
                       scale = FALSE, 
                       na.action = na.fail, 
                       subset = NULL
    )
    # Forward AND backward selection of variables.
    # Selects most significant based on ANOVA and AIC scoring.
    mod.result <- vegan::ordistep(object = mod0, 
                                  scope = stats::formula(mod1),
                                  direction = "both",
                                  trace = FALSE)
    
    # Creates rda model with all significant variables
    mod.signif <- vegan::rda(mod.result$call$formula, data = meta_table)
    
    # Collects most significant columns based on data
    mod.colnames <- stringr::str_match_all(
      string = paste(deparse(mod.result$call$formula[[3]]), collapse = " "),
      pattern = "\\w+" # Takes all words, except "+"
    )
  } else {
    # extracts single column
    mod.colnames = colnames(meta_tab)[grep("RANKSTAT_", names(meta_tab))]
    mod.signif = NULL
  }
  
  result <- list(
    model = mod.signif,
    model.col = mod.colnames
  )
  
  return(result)
}