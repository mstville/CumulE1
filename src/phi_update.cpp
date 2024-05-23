#include "RcppArmadillo.h"
#include "CWVS.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_update(double phi,
                     arma::vec eta,
                     Rcpp::List temp_corr_info,
                     double a_phi,
                     double b_phi,
                     double metrop_var_phi_trans,
                     int acctot_phi_trans){
  int m = eta.size();
  
  arma::mat corr_inv = temp_corr_info[0];
  double log_deter = temp_corr_info[1];
  double psi = log(phi);
  
  double denom = -0.50*log_deter - 0.50*dot(eta, (corr_inv*eta)) + a_phi*psi - b_phi*exp(psi);
  
  double psi_new = R::rnorm(psi, sqrt(metrop_var_phi_trans));
  double phi_new = exp(psi_new);
  Rcpp::List temp_corr_info_new = temporal_corr_fun(m, phi_new);
  arma::mat corr_inv_new = temp_corr_info_new[0];
  double log_deter_new =temp_corr_info_new[1];
  
  double num = -0.50*log_deter_new - 0.50*dot(eta, (corr_inv_new*eta)) + a_phi*psi_new - b_phi*exp(psi_new);
  
  
  double ratio = exp(num - denom);
  int acc = 0;
  if (ratio >= R::runif(0.00, 1.00)){
    phi = phi_new;
    temp_corr_info = temp_corr_info_new;
    acc = 1;
  }
  acctot_phi_trans = acctot_phi_trans + acc;
  
  return Rcpp::List::create(Rcpp::Named("phi") = phi,
                            Rcpp::Named("acctot_phi_trans") = acctot_phi_trans,
                            Rcpp::Named("temp_corr_info") = temp_corr_info);
}
  
