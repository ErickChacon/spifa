
#include <RcppArmadillo.h>
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List spmirt_cpp(Rcpp::NumericVector response,
    int nobs, int nitems, int nfactors, int niter,
    arma::mat constrain_L,
    arma::vec c_initial, arma::vec c_prior_mean, arma::vec c_prior_sd,
    arma::mat A_initial, arma::mat A_prior_mean, arma::mat A_prior_sd,
    arma::vec theta_init
    ) {

  Ifa model(response, nobs, nitems, nfactors,
      constrain_L,
      c_initial, c_prior_mean, c_prior_sd,
      A_initial, A_prior_mean, A_prior_sd,
      theta_init);
  Rcpp::List output = model.sample(niter);
  return output;
}


