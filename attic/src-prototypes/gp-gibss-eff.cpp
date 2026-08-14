
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
Rcpp::List probit_gp_chol(Rcpp::NumericVector y, arma::mat dist, arma::vec params,
    int iter, arma::mat Sigma_proposal) {

  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);

  const int n = dist.n_rows;
  const double sigma2_c = 1;
  arma::vec ones_n(n, arma::fill::ones);
  arma::vec zeros_n(n, arma::fill::zeros);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat X = arma::join_rows(ones_n, eye_n);

  arma::vec params_fix = params;

  // Initializing z
  arma::vec z(n, arma::fill::zeros);
  Rcpp::NumericVector::iterator ity = y.begin();
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  arma::mat z_mat(n, iter);

  // Initializing gp variance covariance
  double tau2 = exp(params(0));
  double phi = exp(params(1));
  arma::mat Sigma_gp(n, n, arma::fill::zeros);
  Sigma_gp = tau2 * exp(- dist / phi);

  //
  arma::mat params_mat(2, iter);


  for (int i = 0; i < iter; ++i) {

    double tau2 = exp(params(0));
    // double tau2 = exp(params_fix(0));
    double phi = exp(params(1));

    // Sampling z
    arma::mat Sigma_theta_prior(n + 1, n + 1, arma::fill::zeros);
    Sigma_theta_prior.submat(0,0, 0,0) = sigma2_c;
    Sigma_theta_prior.submat(1,1, n,n) = Sigma_gp;
    // arma::mat Sigma_theta_post = (X.t() * X + Sigma_theta_prior.i()).i();
    arma::mat Sigma_theta_post =
      arma::inv_sympd(X.t() * X + arma::inv_sympd(Sigma_theta_prior));
    arma::mat S_aux = Sigma_theta_post * X.t();
    arma::vec Mean_theta_post = S_aux * z;

    arma::vec h(n);
    arma::vec w(n);
    arma::vec q(n);
    for (int j = 0; j < n; ++j) {
      h(j) = arma::as_scalar(X.row(j) * S_aux.col(j));
      w(j) = h(j) / (1 - h(j));
      q(j) = sqrt(w(j) + 1);
    }

    ity = y.begin();
    for (int j = 0; j < n; ++j) {
      double z_old = vecsub1(z, j);
      double z_mean = arma::as_scalar(X.row(j) * Mean_theta_post);
      z_mean -= w(j) * (z(j) - z_mean);
      z(j) = RcppTN::rtn1(z_mean, q(j), low_thresh[*ity], high_thresh[*ity]);
      Mean_theta_post += (z(j) - z_old) * S_aux.col(j);
      ity++;
    }

    // Sampling correlation parameters
    arma::mat Sigma_proposal_chol = arma::chol(Sigma_proposal, "lower");
    arma::vec params_aux = params + Sigma_proposal_chol * arma::randn<arma::vec>(2);
    double tau2_aux = exp(params_aux(0));
    // double tau2_aux = exp(params_fix(0));
    double phi_aux = exp(params_aux(1));
    // double phi_aux = phi;

    arma::mat Sigma_z = sigma2_c + Sigma_gp + eye_n;
    arma::mat Sigma_gp_aux(n, n, arma::fill::zeros);
    Sigma_gp_aux = tau2_aux * exp(- dist / phi_aux);
    arma::mat Sigma_z_aux = sigma2_c + Sigma_gp_aux + eye_n;
    // improve line above
    double accept = dmvnorm(z, zeros_n, Sigma_z_aux) +
      R::dnorm(params_aux(0), log(1.0), 0.4, true) +
      R::dnorm(params_aux(1), log(0.03), 0.4, true) -
      dmvnorm(z, zeros_n, Sigma_z) -
      R::dnorm(params(0), log(1.0), 0.4, true) -
      R::dnorm(params(1), log(0.03), 0.4, true);
    if (accept > log(R::runif(0,1))) {
      params = params_aux;
      Sigma_gp = Sigma_gp_aux;
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

