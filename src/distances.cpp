#include "distances.h"
#include "branches.h"

/*-------------------------------------
    Bray-Curtis dissimilarity

    Eq: sum(abs(A - B)) / sum(A + B)
---------------------------------------*/

double BrayCurtis::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    double num = arma::accu(arma::abs(A - B));
    double denum = arma::accu(A + B);
    if (denum == 0.0 || num == 0.0) return 0.0;
    return num / denum;
};

/*----------------------------
    Manhattan Dissimilarity

    Eq. sum(abs(A - B))
-----------------------------*/

double Manhattan::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    return arma::accu(arma::abs(A - B));
};

/*----------------------------
    Jaccard Dissimilarity
-----------------------------*/

double Jaccard::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    double num = arma::accu(arma::min(A,B));
    double denum = arma::accu(arma::max(A,B));
    if (num == 0.0 || denum == 0.0) return 1.0;
    return 1.0 - (num / denum);
};

/*----------------------------
    Cosine Dissimilarity
-----------------------------*/

double Cosine::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    double num = arma::as_scalar(A.t() * B);
    double denomA = arma::norm(A);
    double denomB = arma::norm(B);
    if (denomA == 0.0 || denomB == 0.0 || num == 0.0) return 1.0;
    return 1.0 - (num / (denomA * denomB));
};

/*----------------------------
    Jensen-Shannon Divergence
-----------------------------*/

double JSD::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    arma::mat a = arma::conv_to<arma::mat>::from(A.col(0));
    arma::mat b = arma::conv_to<arma::mat>::from(B.col(0));
    arma::mat m = 0.5 * (a + b);

    // Adding small epsilon for zero's
    arma::mat ai = a + eps;
    arma::mat bi = b + eps;
    arma::mat mi = m + eps;

    arma::vec num = 0.5 * (ai % arma::log(ai / mi));
    arma::vec denum = 0.5 * (bi % arma::log(bi / mi));

    double jsd = arma::accu(num + denum);

    if (jsd == 0.0) return 0.0;
    return jsd;
};


/*----------------------------
    Canberra dissimilarity
-----------------------------*/

double Canberra::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    double sum = 0.0;
    int count = 0;
    arma::sp_mat C = A + B;

    for (auto it = C.begin(); it != C.end(); ++it) {
        int i = it.row();
        int j = it.col();
        double a_val = A(i, j);
        double b_val = B(i, j);
        double num = std::abs(a_val - b_val);
        double denum = std::abs(a_val) + std::abs(b_val);
        if (denum > 0.0) {
            sum += num / denum;
            count++;
        }
    }
    if (count == 0) return std::numeric_limits<double>::quiet_NaN();
    if (sum == 0.0 || count == 0) return 0.0;
    return sum / count;
};

/*----------------------------
    Aitchison dissimilarity
-----------------------------*/

double Aitchison::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    arma::mat a = arma::conv_to<arma::mat>::from(A.col(0));
    arma::mat b = arma::conv_to<arma::mat>::from(B.col(0));

    // Adding small epsilon for zero's
    arma::mat ai = a + eps;
    arma::mat bi = b + eps;

    // Compute geometric means
    double ga = std::exp(arma::as_scalar(arma::mean(arma::log(ai))));
    double gb = std::exp(arma::as_scalar(arma::mean(arma::log(bi))));

    // Compute clr transforms
    arma::vec clr_a = arma::log(ai / ga);
    arma::vec clr_b = arma::log(bi / gb);

    return arma::norm(clr_a - clr_b, 2);
};

/*----------------------------
    UniFrac dissimilarity
-----------------------------*/

UniFrac::UniFrac(const arma::umat& edge_, const arma::vec& edge_lengths_, bool weighted_, bool normalized_)
    : edge(edge_), edge_lengths(edge_lengths_), weighted(weighted_), normalized(normalized_) {}

double UniFrac::operator()(const arma::sp_mat& A, const arma::sp_mat& B) const {
    if (weighted) {
        BranchWeights bw(edge);
        arma::mat weights_mat = bw(A, B);
        arma::vec branch_weights_A = weights_mat.col(0) % edge_lengths;
        arma::vec branch_weights_B = weights_mat.col(1) % edge_lengths;

        double num = arma::accu(arma::abs(branch_weights_A - branch_weights_B));
        if (!normalized) return num;

        double denom = arma::accu(branch_weights_A) + arma::accu(branch_weights_B);
        if (num == 0.0 || denom == 0.0) return 0.0;
        return num / denom;
    } else {
        BranchPresence bp(edge);
        auto presence_pair = bp(A, B);
        const std::vector<bool>& presence_A = presence_pair.first;
        const std::vector<bool>& presence_B = presence_pair.second;

        double distinct = 0.0, shared = 0.0;
        for (size_t i = 0; i < edge.n_rows; ++i) {
            double length = edge_lengths(i);
            bool a = presence_A[i];
            bool b = presence_B[i];
            if (a && b) shared += length;
            else if (a || b) distinct += length;
        }

        double denom = distinct + shared;
        if (distinct == 0.0 || denom == 0.0) return 0.0;
        return distinct / denom;
    }
};