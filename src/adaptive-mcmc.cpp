
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

  const int n = mean.n_elem;
  arma::mat eye_n = arma::eye<arma::mat>(n,n);

  arma::mat Sigma_chol = arma::chol(Sigma, "lower");

  // arma::vec params(n, arma::fill::randn);
  arma::vec params = mean / 3;
  arma::mat params_mat(n, iter);


  for (int i = 0; i < iter; ++i) {

    arma::mat Sigma_proposal(n,n);
    if (i >= (2 * n) && R::runif(0,1) < 0.95) {
      Sigma_proposal = pow(2.38, 2) *
        arma::cov(params_mat.cols(0,i-1).t()) / n;
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

    params_mat.col(i) = params;
  }

  return Rcpp::List::create(
      Rcpp::Named("params") = params_mat.t()
      );
}
