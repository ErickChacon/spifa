
#include <RcppArmadillo.h>
#include "arma-mat.h"
#include "pdf.h"
#include "correlation.h"
#include "name-samples.h"
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
Rcpp::List multi_lm(arma::mat Y, arma::mat X, int iter,
    arma::mat proposal_corr_Sigma, arma::mat proposal_sigmas_Sigma) {

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
  double prior_corr_eta = 1.5;
  double prior_sigmas_mean = 1;
  double prior_sigmas_sd = 0.4;

  // Initializing beta
  arma::vec beta(pq, arma::fill::randn);
  arma::mat Beta = vec2mat(beta, p, q);
  arma::mat beta_samples(pq, iter);

  // Initializing sigmas
  arma::vec sigmas_log_mean(q, arma::fill::randn);
  arma::mat sigmas_log_cov = proposal_sigmas_Sigma;
  arma::vec sigmas_log = sigmas_log_mean;
  arma::vec sigmas = exp(sigmas_log);
  arma::mat sigmas_samples(q, iter);
  double sigmas_logscale = 0;

  // Initializing Corr
  arma::vec corr_free_mean(n_corr, arma::fill::randn);
  arma::mat corr_free_cov = proposal_corr_Sigma;
  arma::vec corr_free = corr_free_mean;
  arma::mat Corr_chol = vec2chol_corr(corr_free, q);
  arma::mat Cov_chol = Corr_chol;
  Cov_chol.each_col() %= sigmas;
  arma::mat Cov_chol_inv = arma::inv(trimatl(Cov_chol));
  arma::mat Cov_inv = Cov_chol_inv.t() * Cov_chol_inv;
  arma::mat corr_samples(round((q+1) * q / 2), iter);
  arma::mat corr_samples2(n_corr, iter);
  double corr_free_logscale = 0;

  // Adaptive stepsize
  const double alpha = 0.8;
  const double C = 0.7;
  double gamma;
  double prob_accept_log;
  double prob_accept_log2;

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
    // arma::mat proposal_corr_Sigma_chol = arma::chol(proposal_corr_Sigma);
    arma::mat proposal_corr_Sigma_chol =
      arma::chol(exp(corr_free_logscale) * corr_free_cov, "lower");
    arma::vec corr_free_aux =
      corr_free + proposal_corr_Sigma_chol * arma::randn<arma::vec>(n_corr);
    arma::mat Corr_chol_aux = vec2chol_corr(corr_free_aux, q);
    arma::mat Cov_chol_aux = Corr_chol_aux;
    Cov_chol_aux.each_col() %= sigmas;
    prob_accept_log =
      dmvnorm_chol(Y.t(), Beta.t() * X.t(), Cov_chol_aux, true) +
      dlkj_corr_free(corr_free_aux, q, prior_corr_eta, true) -
      dmvnorm_chol(Y.t(), Beta.t() * X.t(), Cov_chol, true) -
      dlkj_corr_free(corr_free, q, prior_corr_eta, true);
    if (prob_accept_log > log(R::runif(0,1))) {
      corr_free = corr_free_aux;
      Corr_chol = Corr_chol_aux;
      Cov_chol = Cov_chol_aux;
      Cov_chol_inv = arma::inv(trimatl(Cov_chol));
      Cov_inv = Cov_chol_inv.t() * Cov_chol_inv;
    }

    // Sampling sigmas: random walk Metropolis Hastings
    // arma::mat proposal_sigmas_Sigma_chol = arma::chol(proposal_sigmas_Sigma);
    arma::mat proposal_sigmas_Sigma_chol =
      arma::chol(exp(sigmas_logscale) * sigmas_log_cov);
    arma::vec sigmas_log_aux =
      sigmas_log + proposal_sigmas_Sigma_chol * arma::randn<arma::vec>(q);
    arma::vec sigmas_aux = exp(sigmas_log_aux);
    arma::mat Cov_chol_aux2 = Corr_chol;
    Cov_chol_aux2.each_col() %= sigmas_aux;
    prob_accept_log2 =
      dmvnorm_chol(Y.t(), Beta.t() * X.t(), Cov_chol_aux2, true) +
      dmvnorm(sigmas_aux, arma::ones<arma::vec>(q),
          pow(prior_sigmas_sd, 2) * arma::eye<arma::mat>(q,q)) -
      dmvnorm_chol(Y.t(), Beta.t() * X.t(), Cov_chol, true) -
      dmvnorm(sigmas, arma::ones<arma::vec>(q),
          pow(prior_sigmas_sd, 2) * arma::eye<arma::mat>(q,q));
    if (prob_accept_log2 > log(R::runif(0,1))) {
      sigmas_log = sigmas_log_aux;
      sigmas = sigmas_aux;
    }

    // Update scaling parameter
    gamma = C / pow(i+1, alpha);
    sigmas_logscale += gamma * (std::min(exp(prob_accept_log2), 1.0) - 0.234);
    corr_free_logscale += gamma * (std::min(exp(prob_accept_log), 1.0) - 0.234);

    // Update covariance of parameters with stochastic approximation
    arma::vec center_aux = sigmas_log - sigmas_log_mean;
    sigmas_log_mean += gamma * center_aux;
    sigmas_log_cov += gamma * (center_aux * center_aux.t() - sigmas_log_cov);
    // Update covariance of parameters with stochastic approximation
    arma::vec center_aux2 = corr_free - corr_free_mean;
    corr_free_mean += gamma * center_aux2;
    corr_free_cov += gamma * (center_aux2 * center_aux2.t() - corr_free_cov);

    // Save samples
    beta_samples.col(i) = beta;
    corr_samples.col(i) = trimatl2vec(Corr_chol, true);
    corr_samples2.col(i) = trimatl2vec(Corr_chol * Corr_chol.t(), false);
    sigmas_samples.col(i) = sigmas;
  }

  Rcpp::NumericMatrix beta_samples_rcpp = Rcpp::wrap(beta_samples.t());
  Rcpp::colnames(beta_samples_rcpp) = name_samples_mat(p, q, "Beta");

  Rcpp::NumericMatrix corr_samples_rcpp = Rcpp::wrap(corr_samples.t());
  Rcpp::colnames(corr_samples_rcpp) = name_samples_lower(q, q, "Corr_chol");

  Rcpp::NumericMatrix corr_samples2_rcpp = Rcpp::wrap(corr_samples2.t());
  Rcpp::colnames(corr_samples2_rcpp) = name_samples_lower(q, q, "Corr", false);

  Rcpp::NumericMatrix sigmas_samples_rcpp = Rcpp::wrap(sigmas_samples.t());
  Rcpp::colnames(sigmas_samples_rcpp) = name_samples_vec(q, "sigmas");

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("beta") = beta_samples_rcpp,
      Rcpp::Named("corr_chol") = corr_samples_rcpp,
      Rcpp::Named("corr") = corr_samples2_rcpp,
      Rcpp::Named("sigmas") = sigmas_samples_rcpp
      );

  Rcpp::StringVector myclass(2);
  myclass(0) = "spmirt.list";
  myclass(1) = "list";
  output.attr("class") = myclass;

  return output;

}

