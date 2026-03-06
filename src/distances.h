#ifndef DISTANCES_H
#define DISTANCES_H

#include <RcppArmadillo.h>

struct BrayCurtis {
    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

struct Manhattan { 
    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

struct Euclidean { 
    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

struct Jaccard { 
    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

struct Cosine {
    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

struct JSD { 
    const double eps = 1e-15;
    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

struct Canberra { 
    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

struct UniFrac {
    const arma::umat& edge;
    const arma::vec& edge_lengths;
    bool weighted;
    bool normalized;

    UniFrac(const arma::umat& edge_, const arma::vec& edge_lengths_, bool weighted_, bool normalized_);

    double operator()(const arma::sp_mat& A, const arma::sp_mat& B) const;
};

#endif // DISTANCES_H