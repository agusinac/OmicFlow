unpaired_fold <- function(dt, sample.id, condition_A, condition_B, condition_labels, feature_rank) {
  # subset feature labels before removing them
  feature_labels <- dt[[ feature_rank ]]
  dt <- dt[, .SD, .SDcols = !c(feature_rank)]
  
  
  unpaired_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(unpaired_dt) <- feature_rank
  
  for (i in seq_along(condition_A)) {
    # Subset by condition_A value
    dt_A <- dt[, .SD, .SDcols = colnames(dt)[grepl(condition_A[1], condition.labels)]]
    dt_B <- dt[, .SD, .SDcols = colnames(dt)[grepl(condition_B[1], condition.labels)]]
    
    # Perform subtraction
    dt_diff <- dt_A - dt_B
    
    # Melt into a single column
    dt_diff[, (feature_rank) := feature_labels]
    unpaired_dt <- base::merge(unpaired_dt,
                               data.table::melt(dt_diff,
                                                measure.vars = colnames(dt_diff)[!grepl(feature_rank, colnames(dt_diff))],
                                                variable.name = sample.id, 
                                                value.name = paste0("diff_", i)),
                               by = feature_rank,
                               all.x = TRUE)
    
    
  }
  return(unpaired_dt)
}
