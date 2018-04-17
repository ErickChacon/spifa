
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
#include "pdf.h"
#include "links.h"
// [[Rcpp::depends(RcppTN)]]
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

  const int n = mean.n_elem;
  arma::mat eye_n = arma::eye<arma::mat>(n,n);

  arma::mat Sigma_chol = arma::chol(Sigma, "lower");

  // arma::vec params(n, arma::fill::randn);
  arma::vec params = mean / 3;
  arma::mat params_mat(n, iter);

  arma::vec params_mean(n, arma::fill::zeros);
  arma::mat params_cov(n,n, arma::fill::zeros);


  for (int i = 0; i < iter; ++i) {

    arma::mat Sigma_proposal(n,n);
    if (i >= (2 * n) && R::runif(0,1) < 0.95) {
      Sigma_proposal = pow(2.38, 2) * params_cov / n;
    } else {
      Sigma_proposal = pow(0.1, 2) * arma::eye<arma::mat>(n,n) / n;
    }

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

    if (i > 0) {
      arma::vec params_center = params - params_mean;
      params_cov = params_cov +
        (params_center * params_center.t() - (i+1)/i * params_cov) / (i + 1);
    }
    params_mean = params_mean + (params - params_mean) / (i + 1);

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

  const int n = mean.n_elem;
  arma::mat eye_n = arma::eye<arma::mat>(n,n);

  arma::mat Sigma_chol = arma::chol(Sigma, "lower");

  // arma::vec params(n, arma::fill::randn);
  arma::vec params = mean / 3;
  arma::mat params_mat(n, iter);

  arma::vec params_mean(n, arma::fill::zeros);
  arma::mat params_cov(n,n, arma::fill::zeros);


  for (int i = 0; i < iter; ++i) {

    arma::mat Sigma_proposal(n,n);
    if (i >= (2 * n) && R::runif(0,1) < 0.95) {
      Sigma_proposal = pow(2.38, 2) * params_cov / n;
    } else {
      Sigma_proposal = pow(0.1, 2) * arma::eye<arma::mat>(n,n) / n;
    }

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

    if (i > 0) {
      arma::vec params_center = params - params_mean;
      params_cov = params_cov +
        (params_center * params_center.t() - (i+1)/i * params_cov) / (i + 1);
    }
    params_mean = params_mean + (params - params_mean) / (i + 1);

    params_mat.col(i) = params;
  }

  return Rcpp::List::create(
      Rcpp::Named("params") = params_mat.t()
      );
}

