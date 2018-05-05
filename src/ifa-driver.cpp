
#include <RcppArmadillo.h>
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List spmirt_cpp(
    Rcpp::NumericVector response, arma::mat predictors, arma::mat distances,
    int nobs, int nitems, int nfactors, int niter, int thin,
    arma::mat constrain_L, arma::mat constrain_T, arma::vec constrain_V_sd,
    arma::mat adap_Sigma, double adap_scale, double adap_C,
    double adap_alpha, double adap_accep_prob,
    arma::vec c_initial, arma::vec c_prior_mean, arma::vec c_prior_sd,
    arma::mat A_initial, arma::mat A_prior_mean, arma::mat A_prior_sd,
    arma::mat R_initial, double R_prior_eta,
    arma::mat B_initial, arma::mat B_prior_mean, arma::mat B_prior_sd,
    arma::vec sigmas_gp_initial, arma::vec sigmas_gp_mean, arma::vec sigmas_gp_sd,
    arma::vec phi_gp_initial, arma::vec phi_gp_mean, arma::vec phi_gp_sd,
    std::string model_type
    ) {

  Ifa model(response, predictors, distances,
      nobs, nitems, nfactors,
      constrain_L, constrain_T, constrain_V_sd,
      adap_Sigma, adap_scale,
      c_initial, A_initial, R_initial, B_initial, sigmas_gp_initial, phi_gp_initial,
      model_type);

  Rcpp::List output = model.sample(
      c_prior_mean, c_prior_sd,
      A_prior_mean, A_prior_sd,
      niter, thin,
      adap_C, adap_alpha, adap_accep_prob, R_prior_eta
      );

  return output;
}


