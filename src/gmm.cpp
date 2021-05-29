#include <Rcpp.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace Rcpp;

double mean(NumericVector x){
  //  Compute the mean of a vector x
  double total = 0;
  for(int i = 0; i<x.size(); i++){
    total += x[i];
  }
  return total/x.size();
}

double stdev(NumericVector x){
  // Compute the standard deviation of a vector x
  double var_val = 0;
  double mean_val = mean(x);

  for(int i=0; i<x.size(); i++){
    var_val += (x[i] - mean_val)*(x[i] - mean_val);
  }

  return std::sqrt(1/(x.size()-1.0)*var_val);
}

double dnorm(double x, double mu, double sigma){
  // Returns the normal density at x with mean mu and standard deviation sd
  return std::exp(-std::pow(x-mu, 2.0)/(2.0*std::pow(sigma, 2.0)))*1.0/(sigma*std::sqrt(2.0*M_PI));
}

// Implement a Gaussian Mixture Model in C++
// [[Rcpp::export]]
List gmm(NumericVector X, int d, Nullable<NumericVector> pi_ = R_NilValue,
         Nullable<NumericVector> mu_ = R_NilValue,
         Nullable<NumericVector> sd_ = R_NilValue, int max_iter = 10000, double tol = 1e-5){
  // A classical Gaussian Mixture Model in C++
  // Reference: The Elements of Statistical Learning

  //Initialize pi, mu, sigma
  NumericVector pi (d);
  NumericVector mu (d);
  NumericVector sd (d);

  const int N = X.size();

  //Set pi
  if(pi_.isNull()){
    //Uniform Probability Mass Function
    double temp = 1.0/ (double) d;
    for(int i=0; i<d; i++) pi[i] = temp;
  }else{
    NumericVector temp(pi_);
    pi = temp;
  }

  //Set mu
  if(mu_.isNull()){
    // Random Values
    NumericVector::iterator it = std::min_element(X.begin(), X.end());
    NumericVector::iterator it2 = std::max_element(X.begin(), X.end());
    double Xmin = *it;
    double Xmax = *it2;
    mu = runif(d, Xmin, Xmax);
    std::sort(mu.begin(), mu.end()); //Sort mu
  }else{
    NumericVector temp(mu_);
    mu = temp;
  }

  //Set sigma
  if(sd_.isNull()){
    // Set to SD(X)
    double temp = stdev(X);
    // Rcout << "Standard Deviation" << temp << std::endl;
    for(int i=0; i<d; i++) sd[i] = temp;
  }else{
    NumericVector temp(sd_);
    sd = temp;
  }

  //End Initialization
  //Begin EM algorithm

  // initialize log likelihood
  double old_ll = 0; //for storing the previous iteration log likelihood
  double delta_ll = 1;
  double new_ll = -1e5; //for storing the log likelihood

  NumericMatrix gamma(N,d); //For storing responsibilities

  for(int i = 0; i<max_iter; i++){
    // E-Step
    for(int row = 0; row<N; row++){
      // Loop Through the X_i
      double denominator = 0;
      for(int col=0; col<d; col++){
        // Loop through the pi_i
        gamma(row,col) = pi[col]*dnorm(X[row], mu[col], sd[col]);
        denominator += gamma(row,col);
      }
      for(int col=0; col<d; col++){
        gamma(row,col) /= denominator; //Normalize the gammas
      }
    }

    // M-Step

    for(int col=0; col<d; col++){
      // Compute the mus
      // Reverse the order of these loops!
      mu[col] = 0.0;
      double temp = 0.0;
      for(int row=0; row<N; row++){
        mu[col] += gamma(row,col)*X[row];
        temp += gamma(row,col);
      }
      mu[col] /= temp;
    }

    for(int col=0; col<d; col++){
      //Compute the sigmas
      //Reverse the order of these loops!
      sd[col] = 0.0;
      double temp = 0.0;
      for(int row=0; row<N; row++){
        sd[col] += gamma(row,col)*std::pow(X[row]-mu[col], 2.0);
        temp += gamma(row,col);
      }
      sd[col] /= temp;
      sd[col] = std::sqrt(sd[col]);
    }

    for(int col=0; col<d; col++){
      //Compute the mixing probabilities
      //Reverse the order of these loops!
      pi[col] = 0.0;
      for(int row=0; row<N; row++){
        pi[col] += gamma(row,col)/(double) N;
      }
    }

    // Compute the log-likelihood
    double temp = 0.0;
    new_ll = 0.0;
    for(int row = 0; row<N; row++){
      // Loops in correct order!
      temp = 0.0;
      for(int col = 0; col<d; col++){
        temp += pi[col]*dnorm(X[row], mu[col], sd[col]);
      }
      new_ll += std::log(temp);
    }

    delta_ll = std::abs(new_ll-old_ll);

    // Check Convergence Criteria!
    if(delta_ll < tol && i>=30){
      Rcout << "Converged in " << i << " iterations!" << std::endl;
      return List::create(Named("pi") = pi, Named("mu") =  mu, Named("sd") = sd,
                                Named("loglik") = new_ll);
    }
    old_ll = new_ll;

  }

  Rcout << "Halted after " <<max_iter << " iterations." << std::endl;

  return List::create(Named("pi") = pi, Named("mu") =  mu, Named("sd") = sd,
                      Named("loglik") = new_ll);
}


