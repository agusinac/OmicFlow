#include <RcppArmadillo.h>
#include "distances.h"
#include "PairwiseWorker.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat bray(const arma::sp_mat& mat) {
    BrayCurtis dist;
    return pairwise_dist(mat, dist);
}

// [[Rcpp::export]]
arma::mat jsd(const arma::sp_mat& mat) {
    JSD dist;
    return pairwise_dist(mat, dist);
}

// [[Rcpp::export]]
arma::mat cosine(const arma::sp_mat& mat) {
    Cosine dist;
    return pairwise_dist(mat, dist);
}

// [[Rcpp::export]]
arma::mat jaccard(const arma::sp_mat& mat) {
    Jaccard dist;
    return pairwise_dist(mat, dist);
}

// [[Rcpp::export]]
arma::mat manhattan(const arma::sp_mat& mat) {
    Manhattan dist;
    return pairwise_dist(mat, dist);
}

// [[Rcpp::export]]
arma::mat canberra(const arma::sp_mat& mat) {
    Canberra dist;
    return pairwise_dist(mat, dist);
}

// [[Rcpp::export]]
arma::mat unifrac(const arma::sp_mat& mat, const arma::umat& edge, const arma::vec& edge_lengths, bool weighted, bool normalized) {
    UniFrac dist(edge, edge_lengths, weighted, normalized);
    return pairwise_dist(mat, dist);
};