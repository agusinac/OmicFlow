#include <RcppArmadillo.h>
#include "metrics.h"
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

//------------------------------------------------------------------------//
// Dissimilarity metrics
//------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat bray(const arma::sp_mat& mat) {
    return pairwise_dist(mat, BrayCurtis());
}

// [[Rcpp::export]]
arma::mat euclidean(const arma::sp_mat& mat) {
    return pairwise_dist(mat, Euclidean());
}

// [[Rcpp::export]]
arma::mat cosine(const arma::sp_mat& mat) {
    return pairwise_dist(mat, Cosine());
}

// [[Rcpp::export]]
arma::mat jaccard(const arma::sp_mat& mat) {
    return pairwise_dist(mat, Jaccard());
}

// [[Rcpp::export]]
arma::mat manhattan(const arma::sp_mat& mat) {
    return pairwise_dist(mat, Manhattan());
}

// [[Rcpp::export]]
arma::mat canberra(const arma::sp_mat& mat) {
    return pairwise_dist(mat, Canberra());
}

// [[Rcpp::export]]
arma::mat jsd(const arma::sp_mat& mat) {
    return pairwise_dist(mat, JSD());
}

// // [[Rcpp::export]]
// arma::mat unifrac(const arma::sp_mat& mat, const arma::umat& edge, const arma::vec& edge_lengths, bool weighted, bool normalized) {
//     UniFrac dist(edge, edge_lengths, weighted, normalized);
//     return pairwise_dist(mat, dist);
// };