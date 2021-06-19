#include <stdlib.h>
#include "gmm.h"
#include <random>

using namespace Rcpp ;


//' Mixture Model Bag of Little Bootstraps with Optional Annealing and Non-Uniform
//' Sampling
//'
//' This function returns a list containing the estimated probabilities, means,
//' standard deviations and log likelihoods of the fitted GMM. Implemented in
//' C++
//'
//' @param Y Numeric data vector
//' @param b Bag of Little Bootstraps sample size
//' @param s Number of sample (of size b) to be taken
//' @param r Number of Monte Carlo iterations for each sample
//' @param d Number of distributions in mixture model
//' @param pr_ Optional vector of sampling probabilities. If null, sampling is
//' uniform
//' @param pi_ Optional vector of prior distribution probabilities
//' @param mu_ Optional vector of distribution means
//' @param sd_ Optional vector of distribution standard deviations
//' @param max_iter Maximum number of iterations for EM algorithm. Default is
//' 10,000
//' @param tol Tolerance for convergence of log likelihood for EM algorithm.
//' @param beta Optional parameter representing the starting \eqn{\beta} for
//' annealing. For the default of 1, standard EM algorithm is run. Only needs to
//' be specified for annealing if schedule_ is left NULL
//' @param c Optional secondary parameter that determines growth rate of
//' \eqn{\beta}. Only needs to be specified for annealing if schedule_ is left
//' NULL
//' @param schedule_ Optional scheduling for annealing, if not set to NULL the
//' schedule given will be used in place of the beta and c parameters
//'
//' @return List containing estimates for mixture probabilities, means, and
//' standard deviations
//'
//'
//'
//' @export
// [[Rcpp::export]]


Rcpp::List blb_mix(Rcpp::NumericVector Y, int b, int s, int r, int d,
                   Rcpp::Nullable<Rcpp::NumericVector> pr_ = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> pi_ = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> mu_ = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> sd_ = R_NilValue,
                   int max_iter = 10000, double tol = 1e-5,
                   double beta = 1.0, double c = 1.1,
                   Rcpp::Nullable<Rcpp::NumericVector> schedule_ = R_NilValue){


  Rcpp::NumericMatrix mu_lower(s,d);
  Rcpp::NumericMatrix mu_upper(s,d);
  Rcpp::NumericMatrix sd_lower(s,d);
  Rcpp::NumericMatrix sd_upper(s,d);
  Rcpp::NumericMatrix pi_lower(s,d);
  Rcpp::NumericMatrix pi_upper(s,d);

  Rcpp::NumericMatrix fit_mu(r,d);
  Rcpp::NumericMatrix fit_sd(r,d);
  Rcpp::NumericMatrix fit_pi(r,d);

  int n = Y.length();
  Rcpp::NumericVector pr(n);
  if(pr_.isNull()){
    pr.fill(1); // sample will standardize prob arg
  }else{
    Rcpp::NumericVector pr(n);
    Rcpp::NumericVector temp(pr_);
    pr = temp;
  }

  Rcpp::NumericMatrix piout(s,d);
  Rcpp::NumericMatrix muout(s,d);
  Rcpp::NumericMatrix sdout(s,d);


  for (int j=0; j < s; ++j) {


   // Rcpp::NumericVector samp = RcppArmadillo::sample(index, b, false, pr);
   //Rcpp::NumericVector samp = sample(int n = n, int size = b, bool replace = false, sugar::probs_t probs = pr);
    Rcpp::IntegerVector index = sample(n, b, false, pr);
    Rcpp::NumericVector samp = Y[index];

    Rcpp::NumericMatrix MCpi(r,d);
    Rcpp::NumericMatrix MCmu(r,d);
    Rcpp::NumericMatrix MCsd(r,d);
    Rcpp::NumericVector MCll(r);



    for(int k=0; k < r; ++k) {

      // this sample call's probability could be inversely proportional to pr
      // in case of leverage sampling.
      Rcpp::IntegerVector reindex = sample(b, n, true, pr);
      Rcpp::NumericVector resamp = samp[reindex];


      Rcpp::List out = gmm(Y, d, pi_, mu_, sd_, max_iter, tol, beta, c, schedule_);

      // this could likely be made faster with high d values, but assigning to
      // (_, k) or .column(k) threw Exporter.h line 31 error

      Rcpp::NumericVector pi = out["pi"];
      Rcpp::NumericVector mu = out["mu"];
      Rcpp::NumericVector sd = out["sd"];

      //double loglik = out["loglik"];
      //MCll[r] = loglik;

      for(int i=0; i<d; ++i) {
        MCpi(r,i) = pi[i];
        MCmu(r,i) = mu[i];
        MCsd(r,i) = sd[i];
      }

    }

    Rcpp::NumericVector mean_pi = Rcpp::colMeans(MCpi);
    Rcpp::NumericVector mean_mu = Rcpp::colMeans(MCmu);
    Rcpp::NumericVector mean_sd = Rcpp::colMeans(MCsd);


    // This also likely could be sped up

    for(int i=0; i<d; ++i) {
      piout(s,i) = mean_pi[i];
      muout(s,i) = mean_mu[i];
      sdout(s,i) = mean_sd[i];
    }

  }

  // manually assign data frame attributes to list?
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("pi") = colMeans(piout),
                     Rcpp::Named("mu") = colMeans(muout),
                     Rcpp::Named("sd") = colMeans(sdout));
  out.attr("class") = "data.frame";

  return out;

}
