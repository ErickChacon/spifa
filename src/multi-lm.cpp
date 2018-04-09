
#include <RcppArmadillo.h>
#include "arma-mat.h"
#include "pdf.h"
#include "links.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Multivariate Linear Model
//'
//' @description
//' \code{ifa_gibbs} description.
//'
//' @details
//' details.
//'
//' @param par.
//'
//' @return return.
//'
//' @author Erick A. Chacon-Montalvan
//'
//' @examples
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::List multi_lm(arma::mat Y, arma::mat X, arma::mat Corr, arma::vec sigmas,
    arma::mat Cov, arma::mat Beta, int iter, arma::mat Sigma_proposal) {

  const int n = Y.n_rows;
  const int q = Y.n_cols;
  const int p = X.n_cols;
  const int pq = round(p * q);

  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat eye_q = arma::eye<arma::mat>(q,q);
  arma::mat eye_p = arma::eye<arma::mat>(p,p);
  arma::mat eye_pq = arma::eye<arma::mat>(pq,pq);

  // Priors
  double prior_beta_sigma2 = 1;

  // Initializing beta
  arma::vec beta(pq, arma::fill::randn);
  arma::mat beta_samples(pq, iter);
  // beta_cov_chol(pq, pq);

  for (int i = 0; i < iter; ++i) {

    // Sampling beta

    arma:: mat Cov_inv = arma::inv_sympd(Cov);
    arma::mat beta_cov_inv =
      arma::kron(Cov_inv, X.t() * X) + eye_pq / pow(prior_beta_sigma2, 2);
    arma::mat beta_cov = arma::inv_sympd(beta_cov_inv);
    arma::mat beta_cov_chol = arma::chol(beta_cov, "lower");
    arma::vec beta_mean = beta_cov * arma::vectorise(X.t() * Y * Cov_inv);
    beta = beta_mean + beta_cov_chol * arma::randn<arma::vec>(pq);

    // Save samples
    beta_samples.col(i) = beta;
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = beta_samples.t()
      );

}

