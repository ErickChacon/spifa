
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef IFA_H
#define IFA_H

class Ifa {
private:

  // Metadata
  const int n, q, m, n_corr;
  // Data
  Rcpp::NumericVector y;
  // Parameters to sample
  arma::vec c;                        // difficulty
  arma::vec a;                        // discrimination
  arma::vec theta;                    // latent abilities
  arma::vec z;                        // augmented variable
  // Restrictions
  const arma::mat L;                  // restrictions on discrimination
  // Priors
  const arma::vec c_prior_mean;
  const arma::vec c_prior_sd;
  const arma::vec a_prior_mean;
  const arma::mat A_prior_sd;
  // Parameters-dependent objects
  arma::mat LA;                       // restricted discrimation
  // Constant objects usefull for sampling
  const arma::vec ones_n;
  const arma::mat eye_q;
  const arma::mat eye_n;
  const arma::mat eye_m;
  const arma::vec zeros_nm;
  const Rcpp::NumericVector low_thresh ;
  const Rcpp::NumericVector high_thresh;

public:

  Ifa(Rcpp::NumericVector response, int nobs, int nitems, int nfactors,
      arma::mat constrain_L,
      arma::vec c_ini, arma::vec c_pr_mean, arma::vec c_pr_sd,
      arma::mat A_ini, arma::mat A_pri_mean, arma::mat A_pri_sd,
      arma::vec theta_init
      );

  void update_theta();
  void update_c();
  void update_a();
  void update_z();
  Rcpp::List sample(int niter);
  arma::vec get();
};


#endif /* IFA_H */
