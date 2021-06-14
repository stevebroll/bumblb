#ifndef GMM_H
#define GMM_H
#include <Rcpp.h>

Rcpp::List gmm(Rcpp::NumericVector Y, int s, int r, int d,
               Rcpp::Nullable<Rcpp::NumericVector> pr_ = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> pi_ = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> mu_ = R_NilValue,
               Rcpp::Nullable<Rcpp::NumericVector> sd_ = R_NilValue,
               int max_iter = 10000, double tol = 1e-5,
               double beta = 1.0, double c = 1.1,
               Rcpp::Nullable<Rcpp::NumericVector> schedule_ = R_NilValue);

#endif
