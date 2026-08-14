
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
#include "pdf.h"
#include "links.h"
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
Rcpp::List probit_gp(Rcpp::NumericVector y, arma::mat dist, arma::vec params,
    int iter, arma::mat Sigma_proposal) {

  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);

  const int n = dist.n_rows;
  const double sigma2_c = 1;
  arma::vec ones_n(n, arma::fill::ones);
  arma::vec zeros_n(n, arma::fill::zeros);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat X = arma::join_rows(ones_n, eye_n);

  // Initializing
  arma::vec z(n, arma::fill::zeros);
  Rcpp::NumericVector::iterator ity = y.begin();
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  arma::mat z_mat(n, iter);
  arma::mat params_mat(2, iter);
  arma::vec params_fix = params;

  for (int i = 0; i < iter; ++i) {

    // double tau2 = exp(params(0));
    double tau2 = exp(params_fix(0));
    double phi = exp(params(1));

    // Simulate the auxiliary variables z
    arma::mat Sigma_theta_prior(n + 1, n + 1, arma::fill::zeros);
    Sigma_theta_prior.submat(0,0, 0,0) = sigma2_c;
    Sigma_theta_prior.submat(1,1, n,n) = tau2 * exp(- dist / phi);
    arma::mat Sigma_theta_post =
      arma::inv_sympd(X.t() * X + arma::inv_sympd(Sigma_theta_prior));
    // arma::mat Sigma_theta_post_chol = arma::chol(Sigma_theta_post);
    arma::mat S_aux = Sigma_theta_post * X.t();
    arma::vec Mean_theta_post = S_aux * z;

    arma::vec h(n);
    arma::vec w(n);
    arma::vec q(n);
    for (int j = 0; j < n; ++j) {
      h.subvec(j,j) = X.row(j) * S_aux.col(j);
      w.subvec(j,j) = vecsub1(h, j) / (1 - vecsub1(h, j));
      q.subvec(j,j) = sqrt(vecsub1(w, j) + 1);
    }

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

    // arma::mat Sigma_proposal(2,2, arma::fill::zeros);
    // Sigma_proposal(0,0) = 0.01;
    // Sigma_proposal(1,1) = 0.01;
    // Sigma_proposal(0,1) = 0.005;
    // Sigma_proposal(1,0) = 0.005;
    arma::mat Sigma_proposal_chol = arma::chol(Sigma_proposal, "lower");
    arma::vec params_aux = params + Sigma_proposal_chol * arma::randn<arma::vec>(2);
    // double tau2_aux = exp(params_aux(0));
    double tau2_aux = exp(params_fix(0));
    double phi_aux = exp(params_aux(1));
    // double phi_aux = phi;

    arma::mat Sigma_z = X * Sigma_theta_prior * X.t() + eye_n;
    arma::mat Sigma_theta_prior_aux(n + 1, n + 1, arma::fill::zeros);
    Sigma_theta_prior_aux.submat(0,0, 0,0) = sigma2_c;
    Sigma_theta_prior_aux.submat(1,1, n,n) = tau2_aux * exp(- dist / phi_aux);
    arma::mat Sigma_z_aux = X * Sigma_theta_prior_aux * X.t() + eye_n;
    double accept = dmvnorm(z, zeros_n, Sigma_z_aux) - dmvnorm(z, zeros_n, Sigma_z);
    if (accept > log(R::runif(0,1))) {
      params = params_aux;
    }

    arma::vec params_trans(2);
    params_trans(0) = tau2;
    params_trans(1) = phi;
    z_mat.col(i) = z;
    params_mat.col(i) = params_trans;
  }

  return Rcpp::List::create(
      Rcpp::Named("z") = z_mat.t(),
      Rcpp::Named("param") = params_mat.t()
      );
}

