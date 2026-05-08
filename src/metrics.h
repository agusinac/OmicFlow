#ifndef METRICS_H
#define METRICS_H

#include <RcppArmadillo.h>
#include <cmath>

struct BrayCurtis {
    double num, den;
    BrayCurtis();
    void operator()(double x, double y);
    double distance() const;
};

struct Manhattan {
    double dist;
    Manhattan();
    void operator()(double x, double y);
    double distance() const;
};

struct Euclidean { 
    double sq_dist;
    Euclidean();
    void operator()(double x, double y);
    double distance() const;
};

struct Jaccard {
    double num, den;
    Jaccard();
    void operator()(double x, double y);
    double distance() const;
};

struct Cosine {
    double dot, normA_sq, normB_sq;
    Cosine();
    void operator()(double x, double y);
    double distance() const;
};

struct JSD {
    double num, den;
    JSD();
    void operator()(double x, double y);
    double distance() const;
};

struct Canberra {
    double sum;
    int count;
    Canberra();
    void operator()(double x, double y);
    double distance() const;
};

#endif // METRICS_H