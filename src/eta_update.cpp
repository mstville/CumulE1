#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec eta_update(arma::vec gamma_star,
                     arma::mat corr_inv){
  
  double m = gamma_star.size();
  
  arma::mat cov_eta = inv_sympd(eye(m, m) + corr_inv); 
  
  arma::vec mean_eta = cov_eta*gamma_star;
  
  arma::mat ind_norms = arma::randn(1, m);
  arma::vec eta = mean_eta + trans(ind_norms*arma::chol(cov_eta));
  
  return(eta);
  
}