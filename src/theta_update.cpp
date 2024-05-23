#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double theta_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta,
                    arma::vec gamma,
                    double sig2_theta,
                    double sig2_e){
  
  
  arma::vec zstar=z*gamma;
  double denom = sig2_theta*dot(zstar,zstar) + sig2_e;
  double num = sig2_e*sig2_theta;
  double var_theta = num/denom;
  double mean_theta = (var_theta/sig2_e)*dot(zstar,(y-x*beta));
  
  double theta = R::rnorm(mean_theta, sqrt(var_theta));
  
  return(theta);
}