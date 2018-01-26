
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Spatial Probit Model
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
Rcpp::List probit_sp(arma::mat dist, double tau2, double phi, int iter) {

  const n = dist.n_rows;
  arma::vec ones_n(n, arma::fill::ones);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat X = arma::join_rows(ones_n, eye_n);
  arma::mat Sigma_theta_prior = tau2 * exp(- dist / phi);
  arma::mat Sigma_theta_post = (X.t() * X + Sigma_theta_prior.i()).i();
  arma::mat Sigma_theta_post_chol = arma::chol(Sigma_theta_post);
  arma::mat S_aux = Sigma_theta_post * X.t();

  arma::vec h(n);
  arma::vec w(n);
  arma::vec q(n);

  for (int i = 0; i < n; ++i) {
    h.subvec(i,i) = X.row(i) * S_aux.col(i);
    w.subvec(i,i) = vecsub(h, i) / (1 - vecsub(h, i));
    q.subvec(i,i) = vecsub(w, i) + 1;
  }

  // Initializing
  arma::vec z(n, arma::fill::zeros);
  // Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itmu_z;
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0);
    // ity++;
  }
  arma::vec Mean_theta_post = S_aux * z;

  for (int i = 0; i < iter; ++i) {
    for (int j = 0; int j < n; ++int j) {
      z_old = vecsub(z, j);
      z_mean = arma::as_scalar(X.row(j) * Mean_theta_post);
      z_mean = z_mean -  vecsub(w, j) * (vecsub(z, j) - z_mean);
      z.subvec(j,j) = RcppTN::rtn1(0.0, 1.0);
      Mean_theta_post += (z.subvec(j,j) - z_old) * S_aux.col(j);
    }

    z_mat.col(i) = Mean_theta_post + Sigma_theta_post_chol * arma::randn<arma::vec>(n);
  }

  return Rcpp::List::create(
      Rcpp::Named("theta") = theta_mat.t(),
      Rcpp::Named("c") = c_mat.t(),
      Rcpp::Named("a") = a_mat.t(),
      Rcpp::Named("z") = z_mat.t()
      );
}

