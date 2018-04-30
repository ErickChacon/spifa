
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
#include "name-samples.h"
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

Ifa::Ifa (Rcpp::NumericVector response, int nobs, int nitems, int nfactors,
    arma::mat L_rest, arma::vec theta_init,
    arma::vec c_ini, arma::vec c_pr_mean, arma::vec c_pr_sd,
    arma::mat A_ini, arma::mat A_pri_mean, arma::mat A_pri_sd):
  y(response), n(nobs), q(nitems), m(nfactors),
  n_corr((m-1) * m / 2),
  ones_n(arma::ones(n)),
  zeros_nm(arma::zeros(n*m)),
  eye_q(arma::eye(q,q)),
  eye_n(arma::eye(n,n)),
  eye_m(arma::eye(m,m)),
  low_thresh(Rcpp::NumericVector::create(R_NegInf, 0)),
  high_thresh(Rcpp::NumericVector::create(0, R_PosInf)),
  L(L_rest),
  c_prior_mean(c_pr_mean),
  c_prior_sd(c_pr_sd),
  a_prior_mean(arma::vectorise(A_pri_mean.t())),
  A_prior_sd(A_pri_sd)
{

  // Initializing c: difficulty parameters
  c = arma::randn(q) * 0.3;
  // Initializing a: discrimation parameters
  a = arma::vectorise(A_ini.t());
  LA = L % vec2matt(a, m, q);
  // Initializing z: auxiliary variables
  z = arma::zeros(y.size());
  Rcpp::NumericVector::iterator ity = y.begin();
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  // Initializing theta: latent abilities
  theta = arma::vec(n*m).fill(NA_REAL);

}

void Ifa::update_theta(){
  arma::mat aux_Sigma_inv_chol = arma::chol(LA.t() * LA + eye_m, "lower");
  arma::mat aux_Sigma_chol = arma::inv(trimatl(aux_Sigma_inv_chol)).t();
  arma::mat aux_Sigma = aux_Sigma_chol * aux_Sigma_chol.t();
  arma::mat residual = vec2mat(z - arma::kron(c, ones_n), n, q);
  arma::vec theta_mu = arma::vectorise(residual * LA * aux_Sigma);
  theta = theta_mu + arma::kron(aux_Sigma_chol, eye_n) * arma::randn(n*m);
}

void Ifa::update_c() {
  arma::mat Theta = vec2mat(theta, n, m);
  arma::mat c_X = arma::kron(eye_q, ones_n);
  arma::vec c_Sigma_diag = 1 / (pow(c_prior_sd, -2) + n);
  arma::vec c_mu = c_X.t() * (z - arma::vectorise(Theta * LA.t()));
  c_mu %= c_Sigma_diag;
  c = c_mu + arma::randn(q) % sqrt(c_Sigma_diag);
}

void Ifa::update_a() {
  arma::mat a_Sigma(q*m, q*m, arma::fill::zeros);
  arma::mat a_Sigma_Lt(q*m, q*m, arma::fill::zeros);
  arma::mat a_Sigma_chol(q*m, q*m, arma::fill::zeros);
  arma::mat Theta = vec2mat(theta, n, m);
  arma::mat cross_Theta = Theta.t() * Theta;
  for (int j = 0; j < q; ++j) {
    arma::mat aux_a_Sigma_inv = cross_Theta;
    aux_a_Sigma_inv.rows(find(L.row(j) == 0)) *= 0;
    aux_a_Sigma_inv.cols(find(L.row(j) == 0)) *= 0;
    aux_a_Sigma_inv.diag() += 1 / square(A_prior_sd.row(j));
    // aux_a_Sigma_inv.diag() += 1;
    arma::mat aux_a_Sigma_inv_chol = arma::chol(aux_a_Sigma_inv, "lower");
    arma::mat aux_a_Sigma_chol = arma::inv(trimatl(aux_a_Sigma_inv_chol)).t();
    a_Sigma_chol.submat(j*m, j*m, arma::size(m, m)) = aux_a_Sigma_chol;
    arma::mat aux_a_Sigma = aux_a_Sigma_chol * aux_a_Sigma_chol.t();
    a_Sigma.submat(j*m, j*m, arma::size(m, m)) = aux_a_Sigma;
    // compute Sigma * L.t()
    aux_a_Sigma.cols(find(L.row(j) == 0)) *= 0;
    a_Sigma_Lt.submat(j*m, j*m, arma::size(m, m)) = aux_a_Sigma;
  }
  arma::vec a_mu = a_Sigma_Lt * arma::kron(eye_q, Theta.t()) *
    (z - arma::kron(c, ones_n)) +
    a_Sigma * (a_prior_mean / arma::square(arma::vectorise(A_prior_sd.t())));
  a = a_mu + a_Sigma_chol * arma::randn<arma::vec>(q*m);
  arma::mat A = vec2matt(a, m, q);
  LA = L % A;
}

void Ifa::update_z() {
  arma::mat Theta = vec2mat(theta, n, m);
  arma::vec z_mu = arma::kron(c, ones_n) + arma::vectorise(Theta * LA.t());
  Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itz_mu = z_mu.begin();
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(*itz_mu, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
    itz_mu++;
  }
}

Rcpp::List Ifa::sample(int niter) {

  arma::mat c_samples(q, niter);
  arma::mat a_samples(q*m, niter);
  arma::mat theta_samples(n*m, niter);
  arma::mat z_samples(q*n, niter);

  for (int i = 0; i < niter; ++i) {
    // Update parameters
    update_theta();
    update_c();
    update_a();
    update_z();
    // Save samples
    theta_samples.col(i) = theta;
    c_samples.col(i) = c;
    a_samples.col(i) = arma::vectorise(LA);
    z_samples.col(i) = z;
  }

  Rcpp::NumericMatrix theta_samples_rcpp = Rcpp::wrap(theta_samples.t());
  Rcpp::colnames(theta_samples_rcpp) = name_samples_mat(n, m, "Theta");
  Rcpp::NumericMatrix c_samples_rcpp = Rcpp::wrap(c_samples.t());
  Rcpp::colnames(c_samples_rcpp) = name_samples_vec(q, "c");
  Rcpp::NumericMatrix a_samples_rcpp = Rcpp::wrap(a_samples.t());
  Rcpp::colnames(a_samples_rcpp) = name_samples_mat(q, m, "A");
  Rcpp::NumericMatrix z_samples_rcpp = Rcpp::wrap(z_samples.t());
  Rcpp::colnames(z_samples_rcpp) = name_samples_mat(n, q, "Z");

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("theta") = theta_samples_rcpp,
      Rcpp::Named("c") = c_samples_rcpp,
      Rcpp::Named("a") = a_samples_rcpp,
      Rcpp::Named("z") = z_samples_rcpp
      );

  Rcpp::StringVector myclass(2);
  myclass(0) = "spmirt.list";
  myclass(1) = "list";
  output.attr("class") = myclass;

  return output;
}

arma::vec Ifa::get() {
  return zeros_nm;
}

