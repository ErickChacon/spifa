
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef IFA_H
#define IFA_H

class Ifa {

private:

  // Model type: eifa, cifa, cifa_pred, spifa, spifa_pred
  const std::string model_type;
  // Metadata
  const int n, q, m, p, n_corr;
  // Data
  Rcpp::NumericVector y;              // response variable
  arma::mat dist;                     // distance matrix
  arma::mat X;                        // predictors design matrix
  // Parameters to sample
  arma::vec z;                        // augmented variable
  arma::vec c;                        // difficulty
  arma::vec a;                        // discrimination
  arma::vec theta;                    // latent abilities
  arma::mat B;                        // multivariate fixed effects (pxm)
  arma::vec corr_free;                // V's correlation parameters on free scale
  arma::vec mgp_sd;                   // standard deviations of Gaussian processes
  arma::vec mgp_phi;                  // scale parameters of Gaussian processes
  // Transformed parameters
  arma::mat Corr_chol;                // cholesky of V's correlation matrix
  arma::mat LA;                       // restricted discrimation
  // Restrictions
  const arma::mat L;                  // restrictions on discrimination
  const arma::mat T;                  // restrictions on Gaussian Processes
  const arma::vec V_sd;
  // Adaptive sampling parameters
  arma::vec params;
  arma::vec params_mean;
  arma::mat params_cov;
  double logscale;
  // Hierarquical Priors
  arma::vec theta_prior_mean;
  arma::mat theta_prior_Sigma_chol_inv;
  arma::mat theta_prior_Sigma_inv;
  // Constant objects usefull for sampling
  const arma::vec ones_n;
  const arma::mat eye_q;
  const arma::mat eye_n;
  const arma::mat eye_m;
  const arma::vec zeros_nm;
  const Rcpp::NumericVector low_thresh ;
  const Rcpp::NumericVector high_thresh;

public:

  Ifa(Rcpp::NumericVector response, arma::mat predictors, arma::mat distances,
      int nobs, int nitems, int nfactors,
      arma::mat constrain_L, arma::mat constrain_T, arma::vec constrain_V_sd,
      arma::mat adap_Sigma, double adap_scale,
      arma::vec c_ini, arma::mat A_ini, arma::mat R_ini,
      arma::mat B_ini, arma::vec sigmas_gp_ini, arma::vec phi_gp_ini,
      std::string mod_type
      );

  void update_theta();
  void update_c(const arma::vec& c_prior_mean, const arma::vec& c_prior_sd);
  void update_a(const arma::vec& a_prior_mean, const arma::mat& A_prior_sd);
  void update_z();
  void update_B(const arma::mat& B_prior_sd);
  void update_cov_params(
      const arma::vec sigmas_gp_mean, const arma::vec sigmas_gp_sd,
      const arma::vec phi_gp_mean, const arma::vec phi_gp_sd,
      const double R_prior_eta,
      const double C, const double alpha, const double target, int index);
  Rcpp::List sample(
      int niter, int thin,
      arma::vec c_prior_mean, arma::vec c_prior_sd,
      arma::mat A_prior_mean, arma::mat A_prior_sd,
      double R_prior_eta,
      arma::mat B_prior_mean, arma::mat B_prior_sd,
      arma::vec sigmas_gp_mean, arma::vec sigmas_gp_sd,
      arma::vec phi_gp_mean, arma::vec phi_gp_sd,
      double C, double alpha, double target
      );
};


#endif /* IFA_H */
