
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef IFA_H
#define IFA_H

class Ifa {
private:

  // Model type: eifa, cifa, cifa_pred, spifa, spifa_pred
  const std::string model_type;
  // Metadata
  const int n, q, m, n_corr;
  // Data
  Rcpp::NumericVector y;
  // Parameters to sample
  arma::vec c;                        // difficulty
  arma::vec a;                        // discrimination
  arma::vec theta;                    // latent abilities
  arma::vec z;                        // augmented variable
  arma::vec corr_free;
  // Restrictions
  const arma::mat L;                  // restrictions on discrimination
  // Priors
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

  Ifa(Rcpp::NumericVector response,
      int nobs, int nitems, int nfactors,
      arma::mat constrain_L,
      arma::vec c_ini,
      arma::mat A_ini,
      arma::mat R_ini,
      arma::vec theta_init
      );

  void update_theta();
  void update_c(const arma::vec& c_prior_mean, const arma::vec& c_prior_sd);
  void update_a(const arma::vec& a_prior_mean, const arma::mat& A_prior_sd);
  void update_z();
  void update_cov_params(const double C, const double alpha, const double target,
      const double R_prior_eta);
  Rcpp::List sample(
      arma::vec c_prior_mean, arma::vec c_prior_sd,
      arma::mat A_prior_mean, arma::mat A_prior_sd,
      int niter
      );
};


#endif /* IFA_H */
