## Sparse implementations of distance metrics such as Bray-Curtis, Jaccard Index, Jensen-Shannon divergence and Unifrac
## Functions can be later replaced by cpp sparse implementations.
## TO BE benchmarked with previous tests as done for diversity and feature_glom
## If it works on the 157.000 x 1968 matrix, consider it a success and replacement for rbiom.

ZeroDivision = function(a, b) {
  if (a == 0 | b == 0) {
    return(0)
  } else return(a / b)
}

helperFun = function(sparse_mat, axis = 2) {
  n_samples <- sparse_mat@Dim[axis]

  return(
    list(
      pairs = t(utils::combn(1:n_samples, 2)),
      res_mat = Matrix::Matrix(0,
                               nrow = n_samples,
                               ncol = n_samples,
                               dimnames = list(sparse_mat@Dimnames[[ axis ]], sparse_mat@Dimnames[[ axis ]]),
                               sparse = TRUE)
    )
  )
}

sparseBray = function(sparse_mat, axis = 2) {
  res <- helperFun(sparse_mat = sparse_mat,
                   axis = axis)

  for (i in 1:nrow(res$pairs)) {
    u_i <- res$pairs[i,1]
    v_i <- res$pairs[i,2]

    u <- sparse_mat[, u_i, drop=FALSE]
    v <- sparse_mat[, v_i, drop=FALSE]

    l1_diff <- Matrix::colSums(base::abs(u - v))
    l1_sum <- Matrix::colSums(base::abs(u + v))

    res$res_mat[u_i, v_i] <- ZeroDivision(l1_diff, l1_sum)
    res$res_mat[v_i, u_i] <- res$res_mat[u_i, v_i]
  }

  return(res$res_mat)
}

## This function may be replaced by the more faster implementation: locstra::jaccardMatrix()
## Con: locstra hasnt been updated since 2022 ..
sparseJaccard = function(sparse_mat, axis = 2) {
  res <- helperFun(sparse_mat = sparse_mat,
                   axis = axis)

  # Replace positive numerics to 1
  sparse_mat@x[sparse_mat@x > 0] <- 1

  for (i in 1:nrow(res$pairs)) {
    u_i <- res$pairs[i,1]
    v_i <- res$pairs[i,2]

    u <- sparse_mat[, u_i, drop=FALSE]
    v <- sparse_mat[, v_i, drop=FALSE]

    sum_uv <- Matrix::colSums(u) + Matrix::colSums(v)
    sum_intersection <- Matrix::crossprod(u, v)[1,1]
    sum_union <- sum_uv - sum_intersection

    res$res_mat[u_i, v_i] <- (sum_uv - 2 * sum_intersection) / sum_union
    res$res_mat[v_i, u_i] <- res$res_mat[u_i, v_i]
  }

  return(res$res_mat)
}
