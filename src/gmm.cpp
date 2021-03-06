#include <Rcpp.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <cmath>

//' GMM with Optional Annealing
//'
//' This function returns a list containing the estimated probabilities, means,
//' standard deviations and log likelihoods of the fitted GMM. Implemented in
//' C++
//'
//' @param Y Numeric data vector
//' @param d Number of distributions in mixture model
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
//' schedule given will be used in place of the beta and c parameters.
//'
//' @return List containing probability, standard deviation, and log likelihood
//' values
//'
//'
//' @export
// [[Rcpp::export]]

Rcpp::List gmm(Rcpp::NumericVector Y, int d,
                   Rcpp::Nullable<Rcpp::NumericVector> pi_ = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> mu_ = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> sd_ = R_NilValue,
                   int max_iter = 10000, double tol = 1e-5,
                   double beta = 1.0, double c = 1.1,
                   Rcpp::Nullable<Rcpp::NumericVector> schedule_ = R_NilValue){
  // A classical Gaussian Mixture Model in C++
  // References:  The Elements of Statistical Learning
  //              Deterministic Annealing EM Algorithm
  //              Convergence of the EM algorithm for Gaussian Mixtures with Unbalanced Mixing Coefficients
  //

  //Some assert statements so the code doesn't break
  assert(Y.size()>=3*d); // Make sure there is enough data to fit the gmm
  if(!pi_.isNull()){
    assert(pi_.length()==d);
  }
  if(!mu_.isNull()){
    assert(mu_.length()==d);
  }
  if(!sd_.isNull()){
    assert(sd_.length()==d);
  }

  //Initialize pi, mu, sigma
  Rcpp::NumericVector pi (d);
  Rcpp::NumericVector mu (d);
  Rcpp::NumericVector sd (d);
  Rcpp::NumericVector schedule;

  // Set the annealing schedule
  if(!schedule_.isNull()){
    Rcpp::NumericVector temp(schedule_);
    schedule = temp;
  }
  else if (beta>0 && c>1 && beta<1){
    int N = std::ceil(-std::log(beta)/std::log(c));
    Rcpp::NumericVector temp (N);
    for(int i=0; i<N; i++){
      temp[i] = beta*std::pow(c, i+1);
      if(temp[i]>1){
        temp[i] = 1.0;
      }
    }
    schedule = temp;
  }
  else if (beta == 1){
    Rcpp::NumericVector temp (1);
    temp[0] = 1.0;
    schedule = temp;
  }
  else{
    throw "Improper scheduling criteria.";
  }

  Rcpp::Rcout << "Annealing Schedule:" << std::endl << schedule << std::endl;

  const int N = Y.size();

  //Set pi
  if(pi_.isNull()){
    //Uniform Probability Mass Function
    double temp = 1.0/ (double) d;
    for(int i=0; i<d; i++) pi[i] = temp;
  }else{
    Rcpp::NumericVector temp(pi_);
    pi = temp;
  }

  //Set mu
  if(mu_.isNull()){
    // Random Values
    Rcpp::NumericVector::iterator it  = std::min_element(Y.begin(), Y.end());
    Rcpp::NumericVector::iterator it2 = std::max_element(Y.begin(), Y.end());
    double Ymin = *it;
    double Ymax = *it2;
    mu = Rcpp::runif(d, Ymin, Ymax);
    std::sort(mu.begin(), mu.end()); //Sort mu
  }else{
    Rcpp::NumericVector temp(mu_);
    mu = temp;
  }

  //Set sigma
  if(sd_.isNull()){
    // Set to SD(Y)
    double temp = Rcpp::sd(Y);
    // Rcout << "Standard Deviation" << temp << std::endl;
    for(int i=0; i<d; i++) sd[i] = temp;
  }else{
    Rcpp::NumericVector temp(sd_);
    sd = temp;
  }

  //End Initialization
  //Begin EM algorithm

  // initialize log likelihood
  double old_ll = 0; //for storing the previous iteration log likelihood
  double delta_ll = 1;
  double new_ll = -1e5; //for storing the log likelihood

  Rcpp::NumericMatrix gamma(N,d); //For storing responsibilities
  // double denoms[d];
  for(int bstep = 0; bstep<schedule.size(); bstep++){
    delta_ll = 1;
    old_ll = 0;
    new_ll = -1e5;
    beta = (double) schedule[bstep];

    // Add a small amount of noise
    mu += Rcpp::runif(d,-0.1,0.1); //should make this a function of the data
    // sd += runif(d,-0.1,0.1);
    //std::sort(mu.begin(), mu.end()); //Sort mu

    for(int step = 0; step<max_iter; step++){

      if (step % 50 == 0)
        Rcpp::checkUserInterrupt();

      // E-Step
      for(int row = 0; row<N; row++){
        // Loop Through the Y_i
        double denominator = 0;
        for(int col=0; col<d; col++){
          // Loop through the pi_i
          gamma(row,col) = std::pow(pi[col]*R::dnorm(Y[row], mu[col], sd[col], FALSE),beta);
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
          mu[col] += gamma(row,col)*Y[row];
          temp += gamma(row,col);
        }
        mu[col] /= temp;
      }

      for(int col=0; col<d; col++){
        //Compute the sigmas
        sd[col] = 0.0;
        double temp = 0.0;
        for(int row=0; row<N; row++){
          sd[col] += gamma(row,col)*std::pow(Y[row]-mu[col], 2.0);
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
          temp += pi[col]*R::dnorm(Y[row], mu[col], sd[col], FALSE);
        }
        new_ll += std::log(temp);
      }

      delta_ll = std::abs(new_ll-old_ll);

      // Check Convergence Criteria!
      if(delta_ll < tol && beta==1){
        Rcpp::Rcout << "Converged!" << std::endl;
        return Rcpp::List::create(Rcpp::Named("pi") = pi, Rcpp::Named("mu") =  mu, Rcpp::Named("sd") = sd,
                                  Rcpp::Named("loglik") = new_ll);
      }
      else if(delta_ll < tol){
        break;
      }

      old_ll = new_ll;

    }
    Rcpp::Rcout << beta << std::endl;

  }

  Rcpp::Rcout << "Halted after " <<max_iter << " iterations." << std::endl;

  return Rcpp::List::create(Rcpp::Named("pi") = pi, Rcpp::Named("mu") =  mu, Rcpp::Named("sd") = sd,
                            Rcpp::Named("loglik") = new_ll);
}


