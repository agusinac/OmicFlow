paired_fold <- function(dt, paired.id, sample.id, condition_A, condition_B, unique.id, feature_rank, condition_labels) {
  # tmp data.table
  tmp_dt <- data.table::copy(dt)
  
  feature_labels <- tmp_dt[[ feature_rank ]]
  paired_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(paired_dt) <- feature_rank
  
  for (id in seq_along(unique.id)) {
    # extract samples, takes only matching pairs
    pair <- colnames(dt)[grepl(unique.id[id], colnames(dt))]
    if (length(pair) == 2) {
      sample_pair <- tmp_dt[, .SD, .SDcols = pair]
    } else {
      sample_pair <- NULL
    }
    
    # Computes 2foldChange for matching pairs
    if (!is.null(sample_pair)) {
      # initialize empty table
      dt_diff <- data.table::data.table()
      
      for (i in seq_along(condition_A)) {
        # Subset column pair by condition A and B
        dt_A <- sample_pair[, .SD, .SDcols = colnames(sample_pair)[grepl(condition_A[i], condition_labels)][1]] 
        dt_B <- sample_pair[, .SD, .SDcols = colnames(sample_pair)[grepl(condition_B[i], condition_labels)][1]]
        
        # subtraction
        dt_diff[, (paste0("diff_", i)) := dt_A - dt_B]
        dt_diff[, (feature_rank) := feature_labels]
      }
      dt_diff[, (sample.id) := unique.id[id]]
      paired_dt <- na.omit(rbind(paired_dt, dt_diff, fill = TRUE))
    }
  }
  # removes tmp data table
  base::gc(tmp_dt)
  
  return(paired_dt)
}
