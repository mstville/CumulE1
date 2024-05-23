#ifndef __CWVS__
#define __CWVS__

double sig2_epsilon_update(double a_sig2_e,
                           double b_sig2_e,
                           arma::vec y,
                           arma::mat x,
                           arma::mat z,
                           arma::vec beta,
                           arma::vec alpha);

arma::vec beta_update(arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec y,
                      arma::vec alpha);

double theta_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta,
                    arma::vec gamma,
                    double sig2_theta,
                    double sig2_e);

arma::vec gamma_update(arma::vec y,
                       arma::mat x,
                       arma::mat z,
                       arma::vec beta,
                       arma::vec w,
                       arma::vec gamma_old,
                       arma::vec alpha,
                       double theta,
                       arma::vec eta);

double norm_rs(double a, double b);

double exp_rs(double a, double b);

double half_norm_rs(double a, double b);

double unif_rs(double a, double b);

double rnorm_trunc (double mu, double sigma, double lower, double upper);

arma::vec gamma_star_update(arma::vec gamma,
                            arma::vec eta);

arma::vec eta_update(arma::vec gamma_star,
                     arma::mat corr_inv);

Rcpp::List temporal_corr_fun(int p_z,
                             double phi);

Rcpp::List phi_update(double phi,
                      arma::vec eta,
                      Rcpp::List temp_corr_info,
                      double a_phi,
                      double b_phi,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans);

#endif // __CWVS__