
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
Rcpp::List probit_gp_chol2(Rcpp::NumericVector y, arma::mat dist, arma::vec params,
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
  arma::mat Sigma_marginal = sigma2_c + Sigma_gp;
  arma::mat Sigma_z =  Sigma_marginal;
  Sigma_z.diag() += 1;
  arma::mat Sigma_z_cholinv = arma::inv(trimatl(arma::chol(Sigma_z, "lower")));

  //
  arma::mat params_mat(2, iter);


  for (int i = 0; i < iter; ++i) {

    // Sampling z
    arma::mat Linv_Sigma_marginal = Sigma_z_cholinv * Sigma_marginal;

    arma::vec h(n);
    arma::vec w(n);
    arma::vec q(n);
    for (int j = 0; j < n; ++j) {
      h(j) = Sigma_marginal(j,j) - arma::accu(square(Linv_Sigma_marginal.col(j)));
      w(j) = h(j) / (1 - h(j));
      q(j) = sqrt(w(j) + 1);
    }

    arma::mat Linv_Sigma_marginal2 = Linv_Sigma_marginal.t() * Linv_Sigma_marginal;
    ity = y.begin();
    for (int j = 0; j < n; ++j) {
      double z_mean = arma::as_scalar(Sigma_marginal.row(j) * z);
      z_mean -=
        arma::as_scalar(Linv_Sigma_marginal2.row(j) * z);
      z_mean -= w(j) * (z(j) - z_mean);
      z(j) = RcppTN::rtn1(z_mean, q(j), low_thresh[*ity], high_thresh[*ity]);
      ity++;
    }

    // Sampling correlation parameters
    arma::mat Sigma_proposal_chol = arma::chol(Sigma_proposal, "lower");
    arma::vec params_aux = params + Sigma_proposal_chol * arma::randn<arma::vec>(2);
    // double tau2_aux = exp(params_aux(0));
    double tau2_aux = exp(params_aux(0));
    double phi_aux = exp(params_aux(1));
    // double phi_aux = phi;

    arma::mat Sigma_gp_aux(n, n, arma::fill::zeros);
    Sigma_gp_aux = tau2_aux * exp(- dist / phi_aux);
    arma::mat Sigma_marginal_aux = sigma2_c + Sigma_gp_aux;
    arma::mat Sigma_z_aux = Sigma_marginal_aux + eye_n;
    arma::mat Sigma_z_aux_cholinv =
      arma::inv(trimatl(arma::chol(Sigma_z_aux, "lower")));
    // improve line above
    double accept = dmvnorm_cholinv(z, zeros_n, Sigma_z_aux_cholinv) +
      R::dnorm(params_aux(0), log(1), 0.4, true) +
      R::dnorm(params_aux(1), log(0.03), 0.4, true) -
      dmvnorm_cholinv(z, zeros_n, Sigma_z_cholinv) -
      R::dnorm(params(0), log(1), 0.4, true) -
      R::dnorm(params(1), log(0.03), 0.4, true);
    if (accept > log(R::runif(0,1))) {
      params = params_aux;
      Sigma_gp = Sigma_gp_aux;
      Sigma_z_cholinv = Sigma_z_aux_cholinv;
      Sigma_marginal = Sigma_marginal_aux;
      tau2 = tau2_aux;
      phi = phi_aux;
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

