
#include <RcppArmadillo.h>
#include "arma-mat.h"
#include "links.h"
#include "pdf.h"
#include "correlation.h"
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
Rcpp::List multi_lm(arma::mat Y, arma::mat X, arma::vec sigmas, int iter,
    arma::mat proposal_corr_Sigma) {

  // Define global contants
  const int n = Y.n_rows;
  const int q = Y.n_cols;
  const int p = X.n_cols;
  const int pq = round(p * q);
  const int n_corr = round((q-1) * q / 2);

  // Define global matrices and vectors
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat eye_q = arma::eye<arma::mat>(q,q);
  arma::mat eye_p = arma::eye<arma::mat>(p,p);
  arma::mat eye_pq = arma::eye<arma::mat>(pq,pq);

  // Hyperparameters of priors
  double prior_beta_sigma2 = 10;

  // Initializing beta
  arma::vec beta(pq, arma::fill::randn);
  arma::mat Beta = vec2mat(beta, p, q);
  arma::mat beta_samples(pq, iter);

  // Initializing Corr
  arma::vec corr_free(n_corr);
  arma::mat Corr_chol = vec2chol_corr(corr_free, q);
  arma::mat Cov_chol = Corr_chol;
  Cov_chol.each_col() %= sigmas;
  arma::mat Cov_chol_inv = arma::inv(trimatl(Cov_chol));
  arma::mat Cov_inv = Cov_chol_inv.t() * Cov_chol_inv;
  arma::mat corr_samples(round((q+1) * q / 2), iter);

  for (int i = 0; i < iter; ++i) {

    // Sampling beta: Gibbs
    arma::mat beta_cov_inv =
      arma::kron(Cov_inv, X.t() * X) + eye_pq / pow(prior_beta_sigma2, 2);
    arma::mat beta_cov = arma::inv_sympd(beta_cov_inv);
    arma::mat beta_cov_chol = arma::chol(beta_cov, "lower");
    arma::vec beta_mean = beta_cov * arma::vectorise(X.t() * Y * Cov_inv);
    beta = beta_mean + beta_cov_chol * arma::randn<arma::vec>(pq);
    Beta = vec2mat(beta, p, q);

    // Sampling Corr: random walk Metropolis Hastings
    arma::mat proposal_corr_Sigma_chol = arma::chol(proposal_corr_Sigma);
    arma::vec corr_free_aux =
      corr_free + proposal_corr_Sigma_chol * arma::randn<arma::vec>(n_corr);
    arma::mat Corr_chol_aux = vec2chol_corr(corr_free_aux, q);
    arma::mat Cov_chol_aux = Corr_chol_aux;
    Cov_chol_aux.each_col() %= sigmas;
    double prob_accept_log =
      dmvnorm_chol(Y.t(), Beta.t() * X.t(), Cov_chol_aux, true) +
      dlkj_corr_free(corr_free_aux, q, 1.5, true) -
      dmvnorm_chol(Y.t(), Beta.t() * X.t(), Cov_chol, true) -
      dlkj_corr_free(corr_free, q, 1.5, true);
    if (prob_accept_log > log(R::runif(0,1))) {
      corr_free = corr_free_aux;
      Corr_chol = Corr_chol_aux;
      Cov_chol = Cov_chol_aux;
      Cov_chol_inv = arma::inv(trimatl(Cov_chol));
      Cov_inv = Cov_chol_inv.t() * Cov_chol_inv;
    }

    // Save samples
    beta_samples.col(i) = beta;
    corr_samples.col(i) = trimatl2vec(Corr_chol, true);
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = beta_samples.t(),
      Rcpp::Named("corr_chol") = corr_samples.t()
      );

}

