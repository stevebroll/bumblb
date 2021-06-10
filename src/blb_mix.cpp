#include <Rcpp.h>
#include <stdlib.h>
#include "gmm.h"

// [[Rcpp::export]]

Rcpp::List blb_mix(Rcpp::NumericVector y, int b, int s, int r, int d){
  Rcpp::NumericMatrix mu_lower(s,d);
  Rcpp::NumericMatrix mu_upper(s,d);
  Rcpp::NumericMatrix sd_lower(s,d);
  Rcpp::NumericMatrix sd_upper(s,d);
  Rcpp::NumericMatrix pi_lower(s,d);
  Rcpp::NumericMatrix pi_upper(s,d);
  //const int n = y.size();

  for (int j=0; j < s; ++j) {

    // Sample Here


    Rcpp::NumericMatrix fit_mu(r,d);
    Rcpp::NumericMatrix fit_sd(r,d);
    Rcpp::NumericMatrix fit_pi(r,d);

    for(int k=0; k < r; ++k) {
      // Call gmm.cpp
    }

    // means and CIs here
  }

  // Return list here


  return(0);
}
