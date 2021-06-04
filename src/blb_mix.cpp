#include <Rcpp.h>
#include <stdlib.h>
#include <gmm.h>
using namespace Rcpp;

// [[Rcpp::export]]

List blb_mix(NumericVector y, int b, int s, int r, int d){
  NumericMatrix mu_lower(s,d);
  NumericMatrix mu_upper(s,d);
  NumericMatrix sd_lower(s,d);
  NumericMatrix sd_upper(s,d);
  NumericMatrix pi_lower(s,d);
  NumericMatrix pi_upper(s,d);
  //const int n = y.size();

  for (int j=0; j < s; ++j) {

    // Sample Here


    NumericMatrix fit_mu(r,d);
    NumericMatrix fit_sd(r,d);
    NumericMatrix fit_pi(r,d);

    for(int k=0; k < r; ++k) {
      // Call gmm.cpp
    }

    // means and CIs here
  }

  // Return list here


  return(0);
}
