
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
    arma::mat Cov, arma::mat Beta, int iter, arma::mat proposal_beta_Sigma) {

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

  // Priors
  double prior_beta_sigma2 = 10;

  // Initializing beta
  arma::vec beta(pq, arma::fill::randn);
  arma::mat beta_samples(pq, iter);

  // Initializing Corr
  arma::vec corr_free(n_corr);
  arma::mat Corr_chol = vec2chol_corr(corr_free, q);

  for (int i = 0; i < iter; ++i) {

    // Sampling beta: Gibbs
    arma::mat Cov_inv = arma::inv_sympd(Cov);
    arma::mat beta_cov_inv =
      arma::kron(Cov_inv, X.t() * X) + eye_pq / pow(prior_beta_sigma2, 2);
    arma::mat beta_cov = arma::inv_sympd(beta_cov_inv);
    arma::mat beta_cov_chol = arma::chol(beta_cov, "lower");
    arma::vec beta_mean = beta_cov * arma::vectorise(X.t() * Y * Cov_inv);
    beta = beta_mean + beta_cov_chol * arma::randn<arma::vec>(pq);

    // Sampling Corr: Metropolis Hastings
    arma::vec proposal_beta_Sigma_chol = arma::chol(proposal_beta_Sigma);
    arma::vec corr_free_aux =
      corr_free + proposal_beta_Sigma_chol * arma::randn<arma::vec>(n_corr);
    arma::mat Corr_chol_aux = vec2chol_corr(corr_free_aux, q);
    double prob_accept_log =

      dmvnorm_cholinv(z, zeros_n, Sigma_z_aux_cholinv) +
      R::dnorm(params_aux(0), log(1), 0.4, true) +
      R::dnorm(params_aux(1), log(0.03), 0.4, true) -
      dmvnorm_cholinv(z, zeros_n, Sigma_z_cholinv) -
      R::dnorm(params(0), log(1), 0.4, true) -
      R::dnorm(params(1), log(0.03), 0.4, true);
    if (accept > log(R::runif(0,1))) {
      params = params_aux;
      Sigma_gp = Sigma_gp_aux;
      Sigma_z_cholinv = Sigma_z_aux_cholinv;
      Sigma_marginal = Sigma_marginal_aux;
      tau2 = tau2_aux;
      phi = phi_aux;
    }

    // Save samples
    beta_samples.col(i) = beta;
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = beta_samples.t()
      );

}

