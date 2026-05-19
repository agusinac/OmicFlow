#ifndef UNIFRAC_H
#define UNIFRAC_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/enumerable_thread_specific.h>
#include <cmath>
#include "branches.h"

/*--------------------------------------------------------------------------
  UniFracWorker
--------------------------------------------------------------------------*/
struct UniFracWorker : public RcppParallel::Worker {

    const BranchAccumulator& ba;
    arma::mat&               dist_mat;
    arma::uword              n;
    arma::uword              n_pairs;
    bool weighted;
    bool normalized;

    UniFracWorker(const BranchAccumulator& ba_,
                  arma::mat&               dist_mat_,
                  bool                     weighted_,
                  bool                     normalized_)
        : ba(ba_),
          dist_mat(dist_mat_),
          n(dist_mat_.n_cols),
          n_pairs(dist_mat_.n_cols * (dist_mat_.n_cols - 1) / 2),
          weighted(weighted_),
          normalized(normalized_) {}

    /*----------------------------------------------------------------------
      upper_triangle_ij — maps a flat index k to (i,j) with i<j given precomputed pairs
    ----------------------------------------------------------------------*/
    inline void upper_triangle_ij(std::size_t k,
                                  arma::uword& i, arma::uword& j) const
    {
        i = static_cast<arma::uword>(
                n - 2 - std::floor(
                    std::sqrt(2.0 * (static_cast<double>(n_pairs) - 1.0 - k)
                              + 0.25) - 0.5));
        j = static_cast<arma::uword>(
                k + i + 1 - n_pairs + (n - i) * (n - i - 1) / 2);
    }

    void operator()(std::size_t begin, std::size_t end) override {
        for (std::size_t idx = begin; idx < end; ++idx) {
            arma::uword i, j;
            upper_triangle_ij(idx, i, j);

            double d = weighted
                ? ba.weighted_distance  (i, j, normalized)
                : ba.unweighted_distance(i, j);

            dist_mat(i, j) = d;
            dist_mat(j, i) = d;
        }
    }
};

/*--------------------------------------------------------------------------
  pairwise_unifrac adapted from pairwise_dist
--------------------------------------------------------------------------*/
inline arma::mat pairwise_unifrac(const arma::sp_mat& mat,
                                  const arma::umat&   edge,
                                  const arma::vec&    edge_lengths,
                                  bool                weighted,
                                  bool                normalized)
{
    arma::uword n       = mat.n_cols;
    arma::uword n_pairs = n * (n - 1) / 2;
    arma::mat   dist_mat(n, n, arma::fill::zeros);

    // Tree traversal happens here — once per sample, O(n_samples * n_branches)
    BranchAccumulator ba(edge, edge_lengths, mat);

    UniFracWorker worker(ba, dist_mat, weighted, normalized);
    RcppParallel::parallelFor(0, n_pairs, worker, /*grain=*/256);

    return dist_mat;
}

#endif // UNIFRAC_H