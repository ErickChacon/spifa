
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
Rcpp::List ifa_gibbs(Rcpp::NumericVector y, int n, int q, int N, int m = 1) {

  arma::vec ones_n(n, arma::fill::ones);
  arma::mat eye_q = arma::eye<arma::mat>(q,q);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat eye_m = arma::eye<arma::mat>(m,m);
  arma::mat eye_q_m = arma::eye<arma::mat>(q-m,q-m);

  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);
  // Rcpp::NumericVector lows = low_thresh[y];
  // Rcpp::NumericVector highs = high_thresh[y];

  // Initializing c, a, z

  arma::vec c(q, arma::fill::randn);
  arma::vec mu_c(q);
  arma::mat x_c = arma::kron(eye_q, ones_n);
  arma::mat c_mat(q, N);

  arma::mat A(q,m,arma::fill::zeros);
  A.diag().ones();
  arma::vec a = arma::vectorise(A.t());
  arma::mat Sigma_a((q-m)*m, (q-m)*m);
  arma::mat Sigma_a_aux_chol(m, m);
  arma::mat Sigma_a_aux_chol_inv(m, m);
  arma::mat Sigma_a_chol((q-m)*m, (q-m)*m);
  arma::vec mu_a((q-m)*m);
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
    for (int j = 0; j < m; ++j) {
      // Diagonal parameters
      double sigma2_aux = 1.0 /
        arma::as_scalar(Theta.col(j).t() * Theta.col(j) + 1.0);

      arma::vec bias_aux(n);
      bias_aux.fill(vecsub1(c, j));
      if (j > 1) {
        bias_aux += Theta.cols(0, j-1) * vecsub(a, j * m, j);
      }

      double mean_aux = sigma2_aux *
        arma::as_scalar(Theta.col(j).t() * (vecsub(z, j * n, n) - bias_aux));

      A.submat(j,j, j,j) = RcppTN::rtn1(mean_aux, sqrt(sigma2_aux), 0, R_PosInf);

      // Lower triangular parameters j < m
      if (j > 0) {
        arma::mat aux_eye = arma::eye<arma::mat>(j,j);
        arma::mat low_a_aux_chol = arma::chol(
            Theta.cols(0, j-1).t() * Theta.cols(0,j-1) + aux_eye, "lower");
        arma::mat low_a_aux_chol_inv = arma::inv(low_a_aux_chol);
        arma::mat low_a_Sigma = low_a_aux_chol_inv.t() * low_a_aux_chol_inv;
        arma::mat low_a_Sigma_chol = low_a_aux_chol_inv.t();
        arma::vec low_a_mu = low_a_Sigma * Theta.cols(0,j-1).t() *
          (vecsub(z, j*n, n) - vecsub1(c, j) - Theta.col(j) * matsub1(A, j,j));
        arma::vec low_a = low_a_mu + low_a_Sigma_chol * arma::randn<arma::vec>(j);
        A.submat(j,0, j,j-1) = low_a.t();
      }
    }

    // Rectangle parameters j >= m
    Sigma_a_aux_chol = arma::chol(Theta.t() * Theta + eye_m, "lower");
    Sigma_a_aux_chol_inv = arma::inv(Sigma_a_aux_chol);
    Sigma_a = arma::kron(eye_q_m, Sigma_a_aux_chol_inv.t() * Sigma_a_aux_chol_inv);
    Sigma_a_chol = arma::kron(eye_q_m, Sigma_a_aux_chol_inv.t());
    mu_a = Sigma_a * arma::kron(eye_q_m, Theta.t()) *
      (z.subvec(m*n, q*n-1) - arma::kron(c.subvec(m,q-1), ones_n));
    arma::vec rect_a = mu_a + Sigma_a_chol * arma::randn<arma::vec>((q-m)*m);
    A.submat(m,0, q-1,m-1) = vec2matt(rect_a, m, q-m);
    a = arma::vectorise(A.t());

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

