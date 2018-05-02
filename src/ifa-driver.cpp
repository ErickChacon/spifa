
#include <RcppArmadillo.h>
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List spmirt_cpp(Rcpp::NumericVector response,
    int nobs, int nitems, int nfactors, int niter,
    arma::mat constrain_L,
    arma::vec c_initial, arma::vec c_prior_mean, arma::vec c_prior_sd,
    arma::mat A_initial, arma::mat A_prior_mean, arma::mat A_prior_sd,
    arma::mat R_initial, double R_prior_eta,
    arma::vec theta_init,
    std::string model_type
    ) {

  Ifa model(response,
      nobs, nitems, nfactors,
      constrain_L,
      c_initial,
      A_initial,
      R_initial,
      theta_init);

  Rcpp::List output = model.sample(
      c_prior_mean, c_prior_sd,
      A_prior_mean, A_prior_sd,
      niter);

  return output;
}


