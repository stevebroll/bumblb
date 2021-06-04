#ifndef GMM_H
#define GMM_H
#include <Rcpp.h>

Rcpp::StringVector gmm(Rcpp::NumericVector X, int d,
                                 Rcpp::Nullable<Rcpp::NumericVector> pi_ = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> mu_ = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> sd_ = R_NilValue,
                                 int max_iter = 10000, double tol = 1e-5);
#endif
