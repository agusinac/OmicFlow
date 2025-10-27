#include "branches.h"

arma::mat BranchWeights::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    int n_tips = A.n_rows;
    int n_branches = edge.n_rows;
    arma::mat weights(n_branches, 2);
    arma::vec branch_weights_A(n_branches, 1, arma::fill::zeros);
    arma::vec branch_weights_B(n_branches, 1, arma::fill::zeros);

    // Assign non-zero entries for A
    for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
        int row_idx = it.row();
        int edge_idx = child_to_edge[row_idx];
        if (edge_idx != 1)
            branch_weights_A(edge_idx) = *it;
    }

    // Assign non-zero entries for B
    for (arma::sp_mat::const_iterator it = B.begin(); it != B.end(); ++it) {
        int row_idx = it.row();
        int edge_idx = child_to_edge[row_idx];
        if (edge_idx != 1)
            branch_weights_B(edge_idx) = *it;
    }

    // Sum abundances from tip to root
    for (int i = n_branches - 1; i >= 0; --i) {
        int parent = edge(i, 0);
        if (parent >= n_tips) {
            int pidx = child_to_edge[parent];
            if (pidx != -1) {
                double& bwA_p = branch_weights_A(pidx);
                double& bwA_i = branch_weights_A(i);

                double& bwB_p = branch_weights_B(pidx);
                double& bwB_i = branch_weights_B(i);

                bwA_p += bwA_i;
                bwB_p += bwB_i;
            }
        }
    }
    weights.col(0) = branch_weights_A;
    weights.col(1) = branch_weights_B;
    return weights;
};

std::pair<std::vector<bool>, std::vector<bool>> BranchPresence::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    int n_branches = edge.n_rows;
    int n_tips = A.n_rows;

    std::vector<bool> presence_A(n_branches, false);
    std::vector<bool> presence_B(n_branches, false);

    // Assign non-zero entries for A
    for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
        int idx = it.row();
        int edge_idx = child_to_edge[idx];
        if (edge_idx != -1)
            presence_A[edge_idx] = true;
    }

    // Assign non-zero entries for B
    for (arma::sp_mat::const_iterator it = B.begin(); it != B.end(); ++it) {
        int idx = it.row();
        int edge_idx = child_to_edge[idx];
        if (edge_idx != -1)
            presence_B[edge_idx] = true;
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