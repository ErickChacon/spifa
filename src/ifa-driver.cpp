
#include <RcppArmadillo.h>
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List spmirt_cpp(Rcpp::NumericVector response,
    int nobs, int nitems, int nfactors,
    arma::mat L_rest, int niter, arma::vec theta_init,
    arma::vec c_initial, arma::vec c_prior_mean, arma::vec c_prior_sd,
    arma::mat A_initial, arma::mat A_prior_mean, arma::mat A_prior_sd) {

  Ifa model(response, nobs, nitems, nfactors, L_rest, theta_init,
      c_initial, c_prior_mean, c_prior_sd,
      A_initial, A_prior_mean, A_prior_sd);
  Rcpp::List output = model.sample(niter);
  return output;
}


