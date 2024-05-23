#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_update(arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec y,
                      arma::vec alpha){
  
  int p_x = x.n_cols;
  int n = w.size();
  arma::mat w_mat(n, p_x);
  for(int j = 0; j < p_x; ++j){
     w_mat.col(j) = w;
     }

  arma::mat x_trans = trans(x);

  arma::mat cov_beta = inv_sympd(x_trans*(w_mat%x) + 
                               (1.00/sigma2_beta)*eye(p_x, p_x));
  
  arma::mat L_beta = chol(cov_beta);
  
  arma::vec mean_beta = cov_beta*(x_trans*(w%(y-z*alpha)));
  
  arma::vec ind_norms(p_x); ind_norms.fill(0.00);
  for (int j=0; j<p_x; ++j){
    double zj = R::rnorm(0,1);
    ind_norms(j) = zj;
  }
  
  arma::vec beta = mean_beta + L_beta*ind_norms;
  
  return(beta);
}