#include "metrics.h"

/*-------------------------------------
    Bray-Curtis dissimilarity

    Eq: sum(abs(A - B)) / sum(A + B)
---------------------------------------*/

BrayCurtis::BrayCurtis() : num(0.0), den(0.0) {}
void BrayCurtis::operator()(double x, double y) {
    num += std::abs(x - y);
    den += x + y;
}
double BrayCurtis::distance() const { return den > 0.0 ? num / den : 0.0; }

/*----------------------------
    Manhattan Dissimilarity

    Eq. sum(abs(A - B))
-----------------------------*/

Manhattan::Manhattan() : dist(0.0) {}
void Manhattan::operator()(double x, double y) { dist += std::abs(x - y); }
double Manhattan::distance() const { return dist; }

/*----------------------------
    Euclidean Dissimilarity
-----------------------------*/

Euclidean::Euclidean() : sq_dist(0.0) {}
void Euclidean::operator()(double x, double y) { sq_dist += (x - y) * (x - y); }
double Euclidean::distance() const { return std::sqrt(sq_dist); }


/*----------------------------
    Jaccard Dissimilarity
-----------------------------*/

Jaccard::Jaccard() : num(0.0), den(0.0) {}
void Jaccard::operator()(double x, double y) {
    if (x > 0.0 && y > 0.0) num += std::min(x, y);
    if (x > 0.0 || y > 0.0) den += std::max(x, y);
}
double Jaccard::distance() const { return den > 0.0 ? 1.0 - (num / den) : 1.0; }

/*----------------------------
    Cosine Dissimilarity
-----------------------------*/

Cosine::Cosine() : dot(0.0), normA_sq(0.0), normB_sq(0.0) {}
void Cosine::operator()(double x, double y) {
    dot += x * y;
    normA_sq += x * x;
    normB_sq += y * y;
}
double Cosine::distance() const {
    double na = std::sqrt(normA_sq);
    double nb = std::sqrt(normB_sq);
    return (na == 0.0 || nb == 0.0) ? 1.0 : 1.0 - (dot / (na * nb));
}

/*----------------------------
    Jensen-Shannon Divergence
-----------------------------*/

JSD::JSD() : num(0.0), den(0.0) {}
void JSD::operator()(double x, double y) {
    double m = 0.5 * (x + y);
    if (m > 0.0) {
        if (x > 0.0) num += 0.5 * x * std::log(x / m);
        if (y > 0.0) den += 0.5 * y * std::log(y / m);
    }
}
double JSD::distance() const { 
    return (num + den) / std::log(2.0); 
}

/*----------------------------
    Canberra dissimilarity
-----------------------------*/

Canberra::Canberra() : sum(0.0), count(0) {}
void Canberra::operator()(double x, double y) {
    double numerator = std::abs(x - y);
    double denominator = std::abs(x + y);
    if (denominator > 0.0) {
        sum += numerator / denominator;
        count++;
    }
}
double Canberra::distance() const { 
    return (sum > 0.0 && count > 0) ? sum / count : 0.0; 
}