#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sig2_epsilon_update(double a_sig2_e,
                           double b_sig2_e,
                           arma::vec y,
                           arma::mat x,
                           arma::mat z,
                           arma::vec beta,
                           arma::vec alpha){
  int n = y.size();
  
  double a_sig2_e_update = a_sig2_e + 0.50*n;
  
  double b_sig2_e_update = b_sig2_e + 0.50*dot((y-x*beta-z*alpha),(y-x*beta-z*alpha));
  
  double sig2_epsilon = 1/R::rgamma(a_sig2_e_update, (1.00/b_sig2_e_update));
  
  return(sig2_epsilon);
}