#ifndef BRANCHES_H
#define BRANCHES_H

#include <RcppArmadillo.h>
#include <vector>

struct BranchWeights {
    const arma::umat& edge;
    std::vector<int> child_to_edge;

    BranchWeights(const arma::umat& edge_): edge(edge_) {
        int max_node = edge.max();
        child_to_edge.resize(max_node + 1, -1);
        for (size_t i = 0; i < edge.n_rows; ++i) {
            int child = edge(i, 1);
            child_to_edge[child] = i;
        }
    }
    arma::mat operator()(const arma::sp_mat& tips_A, const arma::sp_mat& tips_B) const;
};

struct BranchPresence {
    const arma::umat& edge;
    std::vector<int> child_to_edge;

    BranchPresence(const arma::umat& edge_) : edge(edge_) {
        int max_node = edge.max();
        child_to_edge.resize(max_node + 1, -1);
        for (size_t i = 0; i < edge.n_rows; ++i) {
            int child = edge(i, 1);
            child_to_edge[child] = i;
        }
    }
    std::pair<std::vector<bool>, std::vector<bool>> operator()(const arma::sp_mat& tips_A, const arma::sp_mat& tips_B) const;
};

#endif // BRANCHES_H