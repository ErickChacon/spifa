
#include <RcppArmadillo.h>
#include "pdf.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Testing Adative Sampling by Haario
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
Rcpp::List adaptive_haario(arma::vec mean, arma::mat Sigma, int iter) {
  // Examples of Adaptive MCMC: Gareth O. Roberts and Jeffrey S. Rosenthal

  // Constants
  const int n = mean.n_elem;
  const arma::mat Sigma_chol = arma::chol(Sigma, "lower");

  // Initializing parameter variables
  arma::vec params_mean(n, arma::fill::zeros);
  arma::mat params_cov(n,n, arma::fill::zeros);
  arma::vec params(n, arma::fill::randn);
  arma::mat params_mat(n, iter);

  for (int i = 0; i < iter; ++i) {

    // Define current Sigma_proposal
    arma::mat Sigma_proposal(n,n);
    if (i >= (2 * n) && R::runif(0,1) < 0.95) {
      Sigma_proposal = pow(2.38, 2) * params_cov / n;
    } else {
      Sigma_proposal = pow(0.1, 2) * arma::eye<arma::mat>(n,n) / n;
    }

    // Update parameters
    arma::mat Sigma_proposal_chol;
    bool nonsingular = arma::chol(Sigma_proposal_chol, Sigma_proposal, "lower");
    if (nonsingular) {
      arma::vec params_aux = params + Sigma_proposal_chol * arma::randn<arma::vec>(n);

      double accept = dmvnorm_chol(params_aux, mean, Sigma_chol, true) -
        dmvnorm_chol(params, mean, Sigma_chol, true);
      if (accept > log(R::runif(0,1))) {
        params = params_aux;
      }
    }

    // Update covariance of parameters recursively
    if (i > 0) {
      arma::vec params_center = params - params_mean;
      params_cov += (params_center*params_center.t() - (i+1)/i*params_cov) / (i+1);
    }
    params_mean += (params - params_mean) / (i + 1);

    // Store current parameter
    params_mat.col(i) = params;
  }

  return Rcpp::List::create(
      Rcpp::Named("params") = params_mat.t()
      );
}

//' @export
// [[Rcpp::export]]
arma::mat variance(arma::mat X) {
  return arma::cov(X);
}


//' @title Testing Adative Sampling by Haario with vanishing adaptation
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
Rcpp::List adaptive_haario_vanish(arma::vec mean, arma::mat Sigma, int iter) {
  // A tutorial on adaptive MCMC: Christophe Andrieu and Johannes Thoms
  // lower alpha leads to faster convergence
  // lower C lead to lower weight at the beginning

  // Constants
  const int n = mean.n_elem;
  const arma::mat Sigma_chol = arma::chol(Sigma, "lower");

  // Initializing parameter variables
  arma::vec params_mean(n, arma::fill::zeros);
  arma::mat params_cov = 2 * arma::eye<arma::mat>(n,n);
  arma::vec params(n, arma::fill::randn);
  arma::mat params_mat(n, iter);

  // Adaptive stepsize
  const double alpha = 0.8; // for 25 parameters
  const double C = 0.7;
  // const double alpha = 0.9; // for 35 parameters
  // const double C = 0.9;
  double gamma;

  for (int i = 0; i < iter; ++i) {

    // Define current Sigma_proposal
    arma::mat Sigma_proposal = pow(2.38, 2) * params_cov / n;

    // Update parameters
    arma::mat Sigma_proposal_chol;
    bool nonsingular = arma::chol(Sigma_proposal_chol, Sigma_proposal, "lower");
    if (nonsingular) {
      arma::vec params_aux = params + Sigma_proposal_chol * arma::randn<arma::vec>(n);

      double accept = dmvnorm_chol(params_aux, mean, Sigma_chol, true) -
        dmvnorm_chol(params, mean, Sigma_chol, true);
      if (accept > log(R::runif(0,1))) {
        params = params_aux;
      }
    }

    // Update covariance of parameters with stochastic approximation
    gamma = C / pow(i+1, alpha);
    arma::vec params_center = params - params_mean;
    params_mean += gamma * params_center;
    params_cov += gamma * (params_center * params_center.t() - params_cov);

    // Store current parameter
    params_mat.col(i) = params;
  }

  return Rcpp::List::create(
      Rcpp::Named("params") = params_mat.t()
      );
}


//' @title AM Algorithm with Global Adaptive Scaling by Christophe
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
Rcpp::List am_vanish_scaling(arma::vec mean, arma::mat Sigma, int iter) {
  // A tutorial on adaptive MCMC: Christophe Andrieu and Johannes Thoms

  // Constants
  const int n = mean.n_elem;
  const arma::mat Sigma_chol = arma::chol(Sigma, "lower");

  // Initializing parameter variables
  arma::vec params_mean(n, arma::fill::zeros);
  arma::mat params_cov = 1 * arma::eye<arma::mat>(n,n);
  arma::vec params(n, arma::fill::randn);
  arma::mat params_mat(n, iter);
  double logscale = 0;

  // Adaptive stepsize
  // As n increases alpha need to go towards 1
  // For n=35, alpha = 0.8 and C = 0.7
  const double alpha = 0.8;
  const double C = 0.7;
  double gamma;
  double accept;

  for (int i = 0; i < iter; ++i) {

    // Define current Sigma_proposal
    arma::mat Sigma_proposal = exp(logscale) * params_cov;

    // Update parameters
    arma::mat Sigma_proposal_chol;
    bool nonsingular = arma::chol(Sigma_proposal_chol, Sigma_proposal, "lower");
    if (nonsingular) {
      arma::vec params_aux = params + Sigma_proposal_chol * arma::randn<arma::vec>(n);

      accept = dmvnorm_chol(params_aux, mean, Sigma_chol, true) -
        dmvnorm_chol(params, mean, Sigma_chol, true);
      if (accept > log(R::runif(0,1))) {
        params = params_aux;
      }
    }

    // Update scaling parameter
    gamma = C / pow(i+1, alpha);
    logscale += gamma * (std::min(exp(accept), 1.0) - 0.234);

    // Update covariance of parameters with stochastic approximation
    arma::vec params_center = params - params_mean;
    params_mean += gamma * params_center;
    params_cov += gamma * (params_center * params_center.t() - params_cov);

    // Store current parameter
    params_mat.col(i) = params;
  }

  return Rcpp::List::create(
      Rcpp::Named("params") = params_mat.t()
      );
}


