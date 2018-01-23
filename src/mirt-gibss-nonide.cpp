
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Bayesian item factor analysis
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
Rcpp::List ifa_gibbs_nonide(Rcpp::NumericVector y, int n, int q, int N, int m = 1) {

  // arma::sp_mat I_q = arma::speye<arma::sp_mat>(5,5);
  // arma::vec ones_n = arma::ones<arma::vec>(n);
  arma::vec ones_n(n, arma::fill::ones);
  arma::mat eye_q = arma::eye<arma::mat>(q,q);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat eye_m = arma::eye<arma::mat>(m,m);
  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);

  // Rcpp::NumericVector lows = low_thresh[y];
  Rcpp::NumericVector lows = low_thresh[y];
  Rcpp::NumericVector highs = high_thresh[y];
  // Rcpp::NumericVector a = rtn(0.0, 1.0, lows, highs);

  // Initializing c, a, z

  arma::vec c(q, arma::fill::randn);
  arma::vec mu_c(q);
  arma::mat x_c = arma::kron(eye_q, ones_n);
  arma::mat c_mat(q, N);

  arma::vec a(q * m, arma::fill::randn);
  arma::mat A = vec2matt(a, m, q);
  arma::mat Sigma_a(q*m, q*m);
  arma::mat Sigma_a_aux_chol(m, m);
  arma::mat Sigma_a_aux_chol_inv(m, m);
  arma::mat Sigma_a_chol(q*m, q*m);
  arma::vec mu_a(q*m);
  arma::mat a_mat(q*m, N);

  arma::vec z = arma::zeros<arma::vec>(y.size());
  arma::vec mu_z(y.size());
  Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itmu_z;
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  arma::mat z_mat(q*n, N);

  arma::vec theta(n*m);
  arma::mat Theta(n, m);
  arma::mat Sigma_theta(n*m, n*m);
  arma::vec mu_theta(n*m);
  arma::mat theta_mat(n*m, N);

  for (int i = 0; i < N; ++i) {

    // Updating latent habilities (theta)
    Sigma_theta = arma::kron(A.t() * A, eye_n) +  arma::eye<arma::mat>(n*m, n*m);
    Sigma_theta = arma::inv(Sigma_theta);
    mu_theta = Sigma_theta * arma::kron(A.t(), eye_n) * (z - arma::kron(c, ones_n));
    theta = mu_theta + arma::chol(Sigma_theta, "lower") * arma::randn<arma::vec>(n*m);
    Theta = vec2mat(theta, n, m);

    // Updating difficulty parameters (c)
    mu_c = x_c.t() * (z - arma::vectorise(Theta * A.t())) / (n + 1);
    c = mu_c + arma::randn<arma::vec>(q) / sqrt(n + 1);

    // Updating discrimation parameters (a)
    Sigma_a_aux_chol = arma::chol(Theta.t() * Theta + eye_m, "lower");
    Sigma_a_aux_chol_inv = arma::inv(Sigma_a_aux_chol);
    Sigma_a = arma::kron(eye_q, Sigma_a_aux_chol_inv.t() * Sigma_a_aux_chol_inv);
    Sigma_a_chol = arma::kron(eye_q, Sigma_a_aux_chol_inv.t());
    mu_a = Sigma_a * arma::kron(eye_q, Theta.t()) * (z - arma::kron(c, ones_n));
    a = mu_a + Sigma_a_chol * arma::randn<arma::vec>(q*m);
    A = vec2matt(a, m, q);

    // Auxiliary variables (z)
    mu_z = x_c * c + arma::kron(eye_q, Theta) * a;
    ity = y.begin();
    itmu_z = mu_z.begin();
    for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
      *it = RcppTN::rtn1(*itmu_z, 1.0, low_thresh[*ity], high_thresh[*ity]);
      ity++;
      itmu_z++;
    }

    // Save samples
    theta_mat.col(i) = theta;
    c_mat.col(i) = c;
    a_mat.col(i) = arma::vectorise(A);
    z_mat.col(i) = z;

  }

  return Rcpp::List::create(
      Rcpp::Named("theta") = theta_mat.t(),
      Rcpp::Named("c") = c_mat.t(),
      Rcpp::Named("a") = a_mat.t(),
      Rcpp::Named("z") = z_mat.t()
      );
}

