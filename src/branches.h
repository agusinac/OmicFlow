#ifndef BRANCHES_H
#define BRANCHES_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

/*--------------------------------------------------------------------------
  BranchAccumulator

  Construction  : O(n_samples * n_branches), tree traversal once per
                  sample
  Per-pair cost : O(n_branches) dot-product / abs-diff loop, no tree walk.

  Weighted normalisation
  --------------------------
  The reference formula is

      d(A,B) = sum_b |w_A[b] - w_B[b]| * len[b]
             / ( sum_b w_A[b]*len[b]  +  sum_b w_B[b]*len[b] )

  where  w_s[b] = (raw subtree abundance of sample s through branch b)
                  / (total abundance of sample s).

  We therefore divide each sample's tip values by its sample total
  *before* propagating up the tree, so that branch_w already stores
  the normalised fractions multiplied by edge length.  The per-pair
  loop then only needs sum |a[b] - b[b]| / (sum a[b] + sum b[b]),
  both of which are trivial dot-products on pre-baked float arrays.

  Tip alignment assumption
  ------------------------
  Rows of the sparse matrix are aligned to tip.labels, meaning row r
  corresponds to the tip whose node index is r.  The edges that have
  child < n_tips therefore have child == row index directly — no extra
  mapping needed.
--------------------------------------------------------------------------*/
struct BranchAccumulator {

    // Tree topology
    std::vector<arma::uword> edge_parent;
    std::vector<arma::uword> edge_child;
    std::vector<double>      edge_lengths_d;
    // Precomputed tree data
    std::vector<int>         parent_edge;
    std::vector<int>         child_to_edge;
    arma::uword              n_branches;

    // Precomputed per-sample data
    std::vector<double> branch_wl;
    std::vector<char>   branch_pres;
    arma::uword         n_samples;

    BranchAccumulator(const arma::umat& edge_,
                      const arma::vec&  edge_lengths_,
                      const arma::sp_mat& mat)
        : n_branches(edge_.n_rows)
    {
        // copy topology into flat arrays
        edge_parent   .resize(n_branches);
        edge_child    .resize(n_branches);
        edge_lengths_d.resize(n_branches);
        for (arma::uword i = 0; i < n_branches; ++i) {
            edge_parent   [i] = edge_(i, 0);
            edge_child    [i] = edge_(i, 1);
            edge_lengths_d[i] = edge_lengths_(i);
        }

        // child_to_edge: node index -> edge row (-1 if not a child)
        arma::uword max_node = *std::max_element(edge_parent.begin(), edge_parent.end());
        max_node = std::max(max_node,
                    *std::max_element(edge_child.begin(), edge_child.end()));
        child_to_edge.assign(max_node + 1, -1);
        for (arma::uword i = 0; i < n_branches; ++i)
            child_to_edge[edge_child[i]] = static_cast<int>(i);

        parent_edge.resize(n_branches);
        for (arma::uword i = 0; i < n_branches; ++i) {
            arma::uword par = edge_parent[i];
            parent_edge[i] = (par < child_to_edge.size())
                             ? child_to_edge[par]
                             : -1;
        }

        precompute(mat);
    }

    /*----------------------------------------------------------------------
      precompute()

      For each sample s:
        1. Read non-zero tips from the sparse column.
        2. Compute sample total (for normalisation).
        3. Seed the per-tip branch weights (normalised).
        4. Propagate up the tree using parent_edge[]
        5. Store branch_wl[s*B .. s*B+B-1] and branch_pres[s*B .. s*B+B-1].
        6. Reset the temporary bw[] buffer via dirty list
    ----------------------------------------------------------------------*/
    void precompute(const arma::sp_mat& mat) {
        n_samples          = mat.n_cols;

        branch_wl  .assign(n_samples * n_branches, 0.0);
        branch_pres.assign(n_samples * n_branches, 0);

        std::vector<double> bw(n_branches, 0.0);
        std::vector<int>    dirty;
        dirty.reserve(n_branches);

        const arma::uword* col_ptrs    = mat.col_ptrs;
        const arma::uword* row_indices = mat.row_indices;
        const double*      values      = mat.values;

        for (arma::uword s = 0; s < n_samples; ++s) {
            arma::uword beg = col_ptrs[s];
            arma::uword fin = col_ptrs[s + 1];

            // sample total for normalisation
            double total = 0.0;
            for (arma::uword k = beg; k < fin; ++k)
                total += values[k];

            double inv_total = (total > 0.0) ? 1.0 / total : 0.0;

            // seed tip edges (tip alignment: row r == node r)
            dirty.clear();
            for (arma::uword k = beg; k < fin; ++k) {
                arma::uword tip  = row_indices[k];
                int         eidx = child_to_edge[tip];
                if (eidx != -1) {
                    bw[eidx] = values[k] * inv_total;
                    dirty.push_back(eidx);
                }
            }

            // propagate up
            for (int i = static_cast<int>(n_branches) - 1; i >= 0; --i) {
                if (bw[i] == 0.0) continue;
                int p = parent_edge[i];
                if (p != -1) {
                    if (bw[p] == 0.0) dirty.push_back(p);
                    bw[p] += bw[i];
                }
            }

            // store into precomputed rows
            double* wl_row   = branch_wl  .data() + s * n_branches;
            char*   pres_row = branch_pres.data() + s * n_branches;
            for (int idx : dirty) {
                wl_row  [idx] = bw[idx] * edge_lengths_d[idx];
                pres_row[idx] = 1;
            }

            // reset only dirty positions
            for (int idx : dirty) bw[idx] = 0.0;
        }
    }

    /*----------------------------------------------------------------------
      weighted_distance()

      d(A,B) = sum_b |wl_A[b] - wl_B[b]|  /  (sum_b wl_A[b] + sum_b wl_B[b])

      wl already has edge length baked in and is normalised by sample total,
      so this exactly matches the reference formula.
      unnormalised case: return the numerator directly.
    ----------------------------------------------------------------------*/
    inline double weighted_distance(arma::uword i, arma::uword j,
                                    bool normalized) const
    {
        const double* __restrict__ a = branch_wl.data() + i * n_branches;
        const double* __restrict__ b = branch_wl.data() + j * n_branches;

        double num = 0.0, sum_a = 0.0, sum_b_val = 0.0;
        for (arma::uword k = 0; k < n_branches; ++k) {
            double ak = a[k], bk = b[k];
            num       += std::abs(ak - bk);
            sum_a     += ak;
            sum_b_val += bk;
        }

        if (!normalized) return num;
        double denom = sum_a + sum_b_val;
        return (denom == 0.0) ? 0.0 : num / denom;
    }

    /*----------------------------------------------------------------------
      unweighted_distance()

      Walks the precomputed presence rows — no tree traversal.
    ----------------------------------------------------------------------*/
    inline double unweighted_distance(arma::uword i, arma::uword j) const
    {
        const char* __restrict__ a = branch_pres.data() + i * n_branches;
        const char* __restrict__ b = branch_pres.data() + j * n_branches;

        double distinct = 0.0, shared = 0.0;
        for (arma::uword k = 0; k < n_branches; ++k) {
            bool ai = a[k], bk = b[k];
            double len = edge_lengths_d[k];
            if      (ai && bk) shared   += len;
            else if (ai || bk) distinct += len;
        }

        double denom = distinct + shared;
        return (denom == 0.0) ? 0.0 : distinct / denom;
    }
};

#endif // BRANCHES_H