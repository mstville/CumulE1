#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec gamma_star_update(arma::vec gamma,
                            arma::vec eta){
  
  arma::vec gamma_star(gamma.size()); gamma_star.fill(0.00);
  int m = gamma.size();
  
  for (int t = 0; t < m; ++t){
    if (gamma(t) == 1.00){
      gamma_star(t) = rnorm_trunc(eta(t), 1.00, 0.00, datum::inf);
    }
    
    if (gamma(t) == 0.00){
      gamma_star(t) = rnorm_trunc(eta(t), 1.00, -datum::inf, 0.00);
    }
  }
  
  return(gamma_star);
}



