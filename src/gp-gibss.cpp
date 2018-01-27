
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
Rcpp::List probit_gp(Rcpp::NumericVector y, arma::mat dist,
    double tau2, double phi, int iter) {

  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);

  const int n = dist.n_rows;
  arma::vec ones_n(n, arma::fill::ones);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat X = arma::join_rows(ones_n, eye_n);
  arma::mat Sigma_theta_prior(n + 1, n + 1, arma::fill::zeros);
  Sigma_theta_prior.submat(0,0, 0,0) = 1;
  Sigma_theta_prior.submat(1,1, n,n) = tau2 * exp(- dist / phi);
  arma::mat Sigma_theta_post = (X.t() * X + Sigma_theta_prior.i()).i();
  arma::mat Sigma_theta_post_chol = arma::chol(Sigma_theta_post);
  arma::mat S_aux = Sigma_theta_post * X.t();

  arma::vec h(n);
  arma::vec w(n);
  arma::vec q(n);

  for (int i = 0; i < n; ++i) {
    h.subvec(i,i) = X.row(i) * S_aux.col(i);
    w.subvec(i,i) = vecsub1(h, i) / (1 - vecsub1(h, i));
    q.subvec(i,i) = sqrt(vecsub1(w, i) + 1);
  }

  // Initializing
  arma::vec z(n, arma::fill::zeros);
  Rcpp::NumericVector::iterator ity = y.begin();
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    // *it = RcppTN::rtn1(50.0, 1.0, R_NegInf, R_PosInf);
    ity++;
  }
  arma::vec Mean_theta_post = S_aux * z;
  arma::mat z_mat(n, iter);

  for (int i = 0; i < iter; ++i) {
    ity = y.begin();
    for (int j = 0; j < n; ++j) {
      double z_old = vecsub1(z, j);
      double z_mean = arma::as_scalar(X.row(j) * Mean_theta_post);
      z_mean -= vecsub1(w, j) * (vecsub1(z, j) - z_mean);
      z.subvec(j,j) =
        RcppTN::rtn1(z_mean, vecsub1(q, j), low_thresh[*ity], high_thresh[*ity]);
        // RcppTN::rtn1(z_mean, vecsub1(q, j), R_NegInf, R_PosInf);
      Mean_theta_post += (vecsub1(z,j) - z_old) * S_aux.col(j);
      ity++;
    }
    z_mat.col(i) = z;
// z_mat.col(i) = Mean_theta_post + Sigma_theta_post_chol * arma::randn<arma::vec>(n);
    // z_mat.col(i) = Mean_theta_post + Sigma_theta_post_chol * arma::randn<arma::vec>(n);
  }

  return Rcpp::List::create(
      Rcpp::Named("z") = z_mat.t()
      // Rcpp::Named("z") = Sigma_theta_post
      );
}

