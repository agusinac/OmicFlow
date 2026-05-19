#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>

/*--------------------------------------------------------------------------
    Template to loop-over-columns and positive indices
----------------------------------------------------------------------------*/

template <typename Metric>
double CSC_iterate(const arma::uword* pos_vec, const arma::uword* row_vec, const double* val_vec,
                     arma::uword sam_i, arma::uword sam_j, Metric metric) {
    
    arma::uword* i = const_cast<arma::uword*>(row_vec + pos_vec[sam_i]);
    arma::uword* j = const_cast<arma::uword*>(row_vec + pos_vec[sam_j]);
    arma::uword* i_end = const_cast<arma::uword*>(row_vec + pos_vec[sam_i + 1]);
    arma::uword* j_end = const_cast<arma::uword*>(row_vec + pos_vec[sam_j + 1]);
    double* val_i = const_cast<double*>(val_vec + pos_vec[sam_i]);
    double* val_j = const_cast<double*>(val_vec + pos_vec[sam_j]);
    
    double x, y;
    
    while (true) {
        if (i != i_end && j != j_end) {
            if (*i == *j) {
                x = *val_i++; i++;
                y = *val_j++; j++;
            } else if (*i < *j) {
                x = *val_i++; i++;
                y = 0.0;
            } else {
                x = 0.0;
                y = *val_j++; j++;
            }
        } else if (i != i_end) {
            x = *val_i++; i++;
            y = 0.0;
        } else if (j != j_end) {
            x = 0.0;
            y = *val_j++; j++;
        } else {
            break;
        }
        metric(x, y);
    }
    return metric.distance();
}

/*----------------------------------------------------------------
    Template for pairwise multithreaded distance computation
----------------------------------------------------------------*/
template <typename Metric>
struct PairwiseDistWorker : public RcppParallel::Worker {
    const arma::uword* pos_vec;
    const arma::uword* row_vec;
    const double* val_vec;
    arma::uword n_cols;
    arma::mat& dist_mat;
    const Metric& metric;

    PairwiseDistWorker(const arma::sp_mat& mat, arma::mat& dist_mat, const Metric& metric)
        :   pos_vec(mat.col_ptrs), 
            row_vec(mat.row_indices),
            val_vec(mat.values), 
            n_cols(mat.n_cols),
            dist_mat(dist_mat), 
            metric(metric) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t idx = begin; idx < end; ++idx) {
            arma::uword sam_i = idx / n_cols;
            arma::uword sam_j = idx % n_cols;
            // only upper triangle
            if (sam_j <= sam_i) continue;
        
            double d = CSC_iterate(pos_vec, row_vec, val_vec, sam_i, sam_j, Metric());

            dist_mat[sam_i * n_cols + sam_j] = d;
        }
    }
};

/*--------------------------------------------------------------------------
    Template for said dissimilarity metric for pairwise samples computation
----------------------------------------------------------------------------*/

template <typename Metric>
arma::mat pairwise_dist(const arma::sp_mat& mat, const Metric& metric) {
    int n = mat.n_cols;
    arma::mat dist_mat(n, n, arma::fill::zeros);

    PairwiseDistWorker<Metric> worker(mat, dist_mat, metric);
    RcppParallel::parallelFor(0, n * n, worker);

    return dist_mat;
};

#endif // UTILS_H