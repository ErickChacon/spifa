
#include <RcppArmadillo.h>
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List spifa_cpp(
    Rcpp::NumericVector response, arma::mat predictors, arma::mat distances,
    int nobs, int nitems, int nfactors, int ngp,
    int niter, int thin, bool standardize,
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
      nobs, nitems, nfactors, ngp,
      constrain_L, constrain_T, constrain_V_sd,
      adap_Sigma, adap_scale,
      c_initial, A_initial, R_initial, B_initial, sigmas_gp_initial, phi_gp_initial,
      model_type);

  Rcpp::List output = model.sample(
      niter, thin, standardize,
      c_prior_mean, c_prior_sd,
      A_prior_mean, A_prior_sd,
      R_prior_eta,
      B_prior_mean, B_prior_sd,
      sigmas_gp_mean, sigmas_gp_sd,
      phi_gp_mean, phi_gp_sd,
      adap_C, adap_alpha, adap_accep_prob
      );

  // Rcpp::List output = Rcpp::List::create(
  //     Rcpp::Named("beta") = 1,
  //     Rcpp::Named("corr_chol") = 2,
  //     Rcpp::Named("corr") = 3,
  //     Rcpp::Named("sigmas") = 5
  //     );
  return output;
}

// [[Rcpp::export]]
Rcpp::List predict_spifa_cpp(arma::mat samples_theta, arma::mat samples_corr_chol,
    arma::mat samples_corr, arma::mat samples_mgp_sd, arma::mat samples_mgp_phi,
    arma::mat samples_betas,
    Rcpp::NumericVector response, arma::mat predictors, arma::mat newpredictors,
    arma::mat distances, arma::mat newdist, arma::mat cross_distances,
    int nobs, int nitems, int nfactors, int ngp, int npred,
    int niter, int burnin, int thin,
    arma::mat constrain_L, arma::mat constrain_T, arma::vec constrain_V_sd,
    std::string model_type
    ) {

  Ifa model(response, predictors, distances,
      nobs, nitems, nfactors, ngp,
      constrain_L, constrain_T, constrain_V_sd,
      model_type);

  Rcpp::List output = model.predict(samples_theta, samples_corr_chol, samples_corr,
      samples_mgp_sd,
      samples_mgp_phi, samples_betas,
      newpredictors, newdist, cross_distances,
      npred, niter, burnin, thin);

  // Rcpp::List output = Rcpp::List::create(
  //     Rcpp::Named("beta") = 1,
  //     Rcpp::Named("corr_chol") = 2,
  //     Rcpp::Named("corr") = 3,
  //     Rcpp::Named("sigmas") = 5
  //     );
  return output;
}

// [[Rcpp::export]]
Rcpp::List predict2_spifa_cpp(arma::mat samples_theta, arma::mat samples_corr_chol,
    arma::mat samples_corr, arma::mat samples_mgp_sd, arma::mat samples_mgp_phi,
    arma::mat samples_betas,
    Rcpp::NumericVector response, arma::mat predictors, arma::mat newpredictors,
    arma::mat distances, arma::mat newdist, arma::mat cross_distances,
    int nobs, int nitems, int nfactors, int ngp, int npred,
    int niter, int burnin, int thin,
    arma::mat constrain_L, arma::mat constrain_T, arma::vec constrain_V_sd,
    std::string model_type
    ) {

  Ifa model(response, predictors, distances,
      nobs, nitems, nfactors, ngp,
      constrain_L, constrain_T, constrain_V_sd,
      model_type);

  Rcpp::List output = model.predict2(samples_theta, samples_corr_chol, samples_corr,
      samples_mgp_sd,
      samples_mgp_phi, samples_betas,
      newpredictors, newdist, cross_distances,
      npred, niter, burnin, thin);

  // Rcpp::List output = Rcpp::List::create(
  //     Rcpp::Named("beta") = 1,
  //     Rcpp::Named("corr_chol") = 2,
  //     Rcpp::Named("corr") = 3,
  //     Rcpp::Named("sigmas") = 5
  //     );
  return output;
}

