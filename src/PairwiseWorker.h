#ifndef PAIRWISEWORKER_H
#define PAIRWISEWORKER_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>

/*----------------------------------------------------------------
    Template for pairwise multithreaded distance computation
----------------------------------------------------------------*/
template <typename Metric>
struct PairwiseDistWorker : public RcppParallel::Worker {
    const arma::sp_mat& mat;
    arma::mat& dist_mat;
    const Metric& dist_fn;

    PairwiseDistWorker(const arma::sp_mat& mat, arma::mat& dist_mat, const Metric& dist_fn)
        : mat(mat), dist_mat(dist_mat), dist_fn(dist_fn) {}

    void operator()(std::size_t begin, std::size_t end) {
        int n = mat.n_cols;
        for (std::size_t idx = begin; idx < end; ++idx) {
            int i = idx / n;
            int j = idx % n;
            if (j <= i) continue; // compute only upper triangle
            arma::sp_mat col_i = mat.col(i);
            arma::sp_mat col_j = mat.col(j);
            double d = dist_fn(col_i, col_j);
            dist_mat(i, j) = d;
            dist_mat(j, i) = d;
        }
    }
};

/*--------------------------------------------------------------------------
    Template for said dissimilarity metric for pairwise samples computation
----------------------------------------------------------------------------*/

template <typename Metric>
arma::mat pairwise_dist(const arma::sp_mat& mat, const Metric& dist_fn) {
    int n = mat.n_cols;
    arma::mat dist_mat(n, n, arma::fill::zeros);

    PairwiseDistWorker<Metric> worker(mat, dist_mat, dist_fn);
    RcppParallel::parallelFor(0, n * n, worker);

    return dist_mat;
}

#endif // PAIRWISEWORKER_H