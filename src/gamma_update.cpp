#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec gamma_update(arma::vec y,
                       arma::mat x,
                       arma::mat z,
                       arma::vec beta,
                       arma::vec w,
                       arma::vec gamma_old,
                       arma::vec alpha,
                       double theta,
                       arma::vec eta){
  
  arma::vec pieces(2); pieces.fill(0.00);
  arma::vec log_pi(2); log_pi.fill(0.00);
  arma::vec gamma = gamma_old;
  arma::vec probs(2);
  int m = z.n_cols;
  
  for (int t = 0; t < m; ++t){
    
    pieces.fill(0.00);
    
    log_pi(0) = log(1.00 - R::pnorm(eta(t),0.00,1.00,true,false));
    
    log_pi(1) = log(R::pnorm(eta(t),0.00,1.00,true,false));
    
    for (int k =0; k <2; ++k){
      gamma(t) = k;
      alpha(t) = gamma(t)*theta;
      pieces(k) = -0.50*dot((y-x*beta-z*alpha), w%(y-x*beta-z*alpha)) + log_pi(k);
    }
    
    probs.fill(0.00);
    
    for (int k = 0; k <2; ++k){
      probs(k) = 1.00/sum(exp(pieces - pieces(k)));
      
      if(arma::is_finite(probs(k)) == 0.00){
        probs(k) = 0.00;  /*Computational Correction*/
      }
    }
    
    gamma(t) = as<double>(Rcpp::rbinom(1, 1, probs(1)));
    alpha(t) = gamma(t)*theta;
    
    if((gamma(t) == 0) && (t < m - 1)){
      for(int l = t+1; l < m; ++l){
        gamma(l) = 0;
      }
      t = m;
    }
  }
  
  return(gamma);
  
}