#include <Rcpp.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace Rcpp;

// Implement a Gaussian Mixture Model in C++
// [[Rcpp::export]]
List gmm_annealing(NumericVector X, int d,
         Nullable<NumericVector> pi_ = R_NilValue,
         Nullable<NumericVector> mu_ = R_NilValue,
         Nullable<NumericVector> sd_ = R_NilValue,
         int max_iter = 10000, double tol = 1e-5,
         double beta = 0.1, double c = 1.1){
  // A classical Gaussian Mixture Model in C++
  // References:  The Elements of Statistical Learning
  //              Deterministic Annealing EM Algorithm
  // To Do: Add the ability to input an annealing schedule

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
    double temp = Rcpp::sd(X);
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
  // double denoms[d];
  beta /=c;
  do{
    delta_ll = 1;
    old_ll = 0;
    new_ll = -1e5;
    beta *= c;
    if(beta>=1) beta=1;

    // Add a small amount of noise
    mu += runif(d,-0.1,0.1); //should make this a function of the data
    // sd += runif(d,-0.1,0.1);
    //std::sort(mu.begin(), mu.end()); //Sort mu

    for(int step = 0; step<max_iter; step++){

      if (step % 50 == 0)
        Rcpp::checkUserInterrupt();

      // E-Step
      for(int row = 0; row<N; row++){
        // Loop Through the X_i
        double denominator = 0;
        for(int col=0; col<d; col++){
          // Loop through the pi_i
          gamma(row,col) = std::pow(pi[col]*R::dnorm(X[row], mu[col], sd[col], FALSE),beta);
          denominator += gamma(row,col);
        }
        for(int col=0; col<d; col++){
          gamma(row,col) /= denominator; //Normalize the gammas
        }
      }

      // M-Step

      for(int col=0; col<d; col++){
        // Compute the mus
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
        sd[col] = 0.0;
        double temp = 0.0;
        for(int row=0; row<N; row++){
          sd[col] += gamma(row,col)*std::pow(X[row]-mu[col], 2.0);
          temp += gamma(row,col);
        }
        sd[col] /= temp;
        sd[col] = std::sqrt(sd[col]);
      }

      //Reset pi before we set it
      for(int col = 0; col<d; col++) pi[col] = 0.0;

      //Compute pi
      for(int row = 0; row<N; row++){
        for(int col = 0;col<d; col++){
          pi[col]+= gamma(row,col)/(double) N;
        }
      }

      // Compute the log-likelihood
      double temp = 0.0;
      new_ll = 0.0;
      for(int row = 0; row<N; row++){
        // Loops in correct order!
        temp = 0.0;
        for(int col = 0; col<d; col++){
          temp += pi[col]*R::dnorm(X[row], mu[col], sd[col], FALSE);
        }
        new_ll += std::log(temp);
      }

      delta_ll = std::abs(new_ll-old_ll);

      // Check Convergence Criteria!
      if(delta_ll < tol && beta==1){
        Rcout << "Converged!" << std::endl;
        return List::create(Named("pi") = pi, Named("mu") =  mu, Named("sd") = sd,
                                  Named("loglik") = new_ll);
      }
      else if(delta_ll < tol){
        break;
      }

      old_ll = new_ll;

    }
    Rcout << beta << std::endl;

  } while (beta<1);

  Rcout << "Halted after " <<max_iter << " iterations." << std::endl;

  return List::create(Named("pi") = pi, Named("mu") =  mu, Named("sd") = sd,
                            Named("loglik") = new_ll);
}


