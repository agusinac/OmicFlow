unpaired_fold <- function(dt, sample.id, condition_A, condition_B, condition_labels, feature_rank, cpus = 8) {
  # Creates tmp data table
  tmp_dt <- data.table::copy(dt)
  
  # subset feature labels before removing them
  feature_labels <- tmp_dt[[ feature_rank ]]
  tmp_dt <- tmp_dt[, .SD, .SDcols = !c(feature_rank)]
  
  # Create data.tables for results
  unpaired_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(unpaired_dt) <- feature_rank
  
  volcano_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(volcano_dt) <- feature_rank
  
  # Register Parallel backend
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  
  # Computing for multiple conditions
  for (i in seq_along(condition_A)) {
    # Subset by condition_A value
    dt_A <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_A[i], condition_labels)]]
    dt_B <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_B[i], condition_labels)]]
    
    # Improve with foreach and parallelize it!
    # Create cross-wise combinations
    combinations <- data.table::data.table(expand.grid(colnames(dt_A), colnames(dt_B)))
    
    # Perform subtraction cross wise
    results <- foreach(j = seq_along(combinations$Var1), .combine = cbind, .packages = 'data.table') %dopar% {
      col_A <- combinations$Var1[j]
      col_B <- combinations$Var2[j]
      
      dt_A[, ..col_A] - dt_B[, ..col_B]
    }
    
    unpaired_dt <- cbind(unpaired_dt, results)
    
    # Compute pvalues with wilcox test
    mat_A <- as.matrix(dt_A)
    mat_B <- as.matrix(dt_B)
    for (k in seq_along(feature_labels)) {
      # save p-values in data.table
      volcano_dt[k, (paste0("pvalue_", i)) := stats::wilcox.test(mat_A[k, ], mat_B[k, ], correct = TRUE)$p.value]
    }
    
    # Compute row means for each taxa
    volcano_dt[, (paste0("foldchange_", i)) := rowMeans(unpaired_dt[, .SD, .SDcols = !c(feature_rank)])]
    
    # Melt into a single column
    final_dt <- data.table::melt(unpaired_dt,
                                 id.vars = feature_rank,
                                 variable.name = sample.id, 
                                 value.name = paste0("diff_", i))
  }
  # Stop the cluster
  parallel::stopCluster(cl)
  
  result <- list(
    data = final_dt,
    volcano = volcano_dt
  )
  return(result)
}
