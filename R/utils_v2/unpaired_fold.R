unpaired_fold <- function(dt, sample.id, condition_A, condition_B, condition_labels, feature_rank) {
  # subset feature labels before removing them
  feature_labels <- dt[[ feature_rank ]]
  dt <- dt[, .SD, .SDcols = !c(feature_rank)]
  
  # Create data.tables for results
  unpaired_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(unpaired_dt) <- feature_rank
  
  pvalues_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(pvalues_dt) <- feature_rank
  
  volcano_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(volcano_dt) <- feature_rank

  # Computing for multiple conditions
  for (i in seq_along(condition_A)) {
    # Subset by condition_A value
    dt_A <- dt[, .SD, .SDcols = colnames(dt)[grepl(condition_A[i], condition_labels)]]
    dt_B <- dt[, .SD, .SDcols = colnames(dt)[grepl(condition_B[i], condition_labels)]]
    
    # Improve with foreach and parallelize it!
    # Create cross-wise combinations
    combinations <- data.table::data.table(expand.grid(colnames(dt_A), colnames(dt_B)))
    
    # Perform subtraction cross wise
    for (j in seq_along(combinations$Var1)) {
      col_A <- combinations$Var1[j]
      col_B <- combinations$Var2[j]
      
      unpaired_dt <- cbind(unpaired_dt, dt_A[, ..col_A] - dt_B[, ..col_B])
    }
    
    # Compute pvalues with wilcox test
    mat_A <- as.matrix(dt_A)
    mat_B <- as.matrix(dt_B)
    for (k in seq_along(feature_labels)) {
      # save p-values in data.table
      pvalues_dt[k, "pvalue" := stats::wilcox.test(mat_A[k, ], mat_B[k, ], correct = TRUE)$p.value]
      
      # save results in volcano data.table
      volcano_dt[k, "pvalue" := pvalues_dt[k, "pvalue"]]
    }
    
    # Compute row means for each taxa
    unpaired_dt[, "Mean" := rowMeans(.SD), by = "Genus"]
    volcano_dt$Mean <- unpaired_dt$Mean
    
    # Melt into a single column
    final_dt <- data.table::melt(unpaired_dt[, .SD, .SDcols = !c("Mean")],
                                 id.vars = feature_rank,
                                 variable.name = sample.id, 
                                 value.name = paste0("diff_", i))
  }
  result <- list(
    data = final_dt,
    pvalues = pvalues_dt,
    volcano = volcano_dt
  )
  
  return(result)
}
