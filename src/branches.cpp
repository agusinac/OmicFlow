#include "branches.h"

arma::mat BranchWeights::operator()(const arma::sp_mat& tips_A, const arma::sp_mat& tips_B) const {
    int n_branches = edge.n_rows;
    int n_tips = tips_A.n_rows;
    arma::mat weights(n_branches, 2);
    arma::vec branch_weights_A(n_branches, 1, arma::fill::zeros);
    arma::vec branch_weights_B(n_branches, 1, arma::fill::zeros);

    // Assign tip abundances to branches using direct lookup
    for (int i = 0; i < n_tips; ++i) {
        int edge_idx = child_to_edge[i];
        if (edge_idx == -1) continue;
        branch_weights_A(edge_idx) = tips_A(i, 0);
        branch_weights_B(edge_idx) = tips_B(i, 0);
    }
    // Propagate abundances upward (same as before)
    // Error `Mat::operator(): Index out of bounds` happens at this code snippet!
    for (int i = n_branches - 1; i >= 0; --i) {
        int parent = edge(i, 0);
        if (parent >= n_tips) {
            int pidx = child_to_edge[parent];
            if (pidx != -1) {
                branch_weights_A(pidx) += branch_weights_A(i);
                branch_weights_B(pidx) += branch_weights_B(i);
            }
        }
    }
    weights.col(0) = branch_weights_A;
    weights.col(1) = branch_weights_B;
    return weights;
};

std::pair<std::vector<bool>, std::vector<bool>> BranchPresence::operator()(const arma::sp_mat& tips_A, const arma::sp_mat& tips_B) const {
    int n_branches = edge.n_rows;
    int n_tips = tips_A.n_rows;

    std::vector<bool> presence_A(n_branches, false);
    std::vector<bool> presence_B(n_branches, false);

    for (int i = 0; i < n_tips; ++i) {
        int edge_idx = child_to_edge[i];
        if (edge_idx != -1) {
            presence_A[edge_idx] = (tips_A(i, 0) > 0);
            presence_B[edge_idx] = (tips_B(i, 0) > 0);
        }
    }

    for (int i = n_branches - 1; i >= 0; --i) {
        int parent = edge(i, 0);
        if (parent >= n_tips) {
            int pidx = child_to_edge[parent];
            if (pidx != -1) {
                presence_A[pidx] = presence_A[pidx] || presence_A[i];
                presence_B[pidx] = presence_B[pidx] || presence_B[i];
            }
        }
    }

    return {presence_A, presence_B};
};