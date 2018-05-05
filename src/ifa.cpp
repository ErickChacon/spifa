
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
#include "name-samples.h"
#include "correlation.h"
#include "ifa.h"
#include "pdf.h"
// [[Rcpp::depends(RcppArmadillo)]]

Ifa::Ifa (Rcpp::NumericVector response, arma::mat predictors, arma::mat distances,
    int nobs, int nitems, int nfactors,
    arma::mat constrain_L, arma::mat constrain_T, arma::vec constrain_V_sd,
    arma::mat adap_Sigma, double adap_scale,
    arma::vec c_ini, arma::mat A_ini, arma::mat R_ini,
    arma::mat B_ini, arma::vec sigmas_gp_ini, arma::vec phi_gp_ini,
    std::string mod_type
    ):
  model_type(mod_type),
  y(response), n(nobs), q(nitems), m(nfactors),
  n_corr((m-1) * m / 2),
  ones_n(arma::ones(n)),
  zeros_nm(arma::zeros(n*m)),
  eye_q(arma::eye(q,q)),
  eye_n(arma::eye(n,n)),
  eye_m(arma::eye(m,m)),
  low_thresh(Rcpp::NumericVector::create(R_NegInf, 0)),
  high_thresh(Rcpp::NumericVector::create(0, R_PosInf)),
  L(constrain_L),
  T(constrain_T),
  V_sd(constrain_V_sd),
  logscale(log(adap_scale)),
  params_cov(adap_Sigma)
{

  // Initializing c: difficulty parameters
  c = c_ini;
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
  // Initializing correlation parameters
  Corr_chol = arma::chol(R_ini, "lower");
  corr_free = chol_corr2vec(Corr_chol);
  arma::mat Cov_chol = Corr_chol;
  Cov_chol.each_col() %= V_sd;

  // Join covariance parameters
  // params = arma::join_cols(mgp_sd_log, mgp_phi_log);
  // params = arma::join_cols(params, corr_free);
  if (model_type != "spifa" && model_type != "spifa_pred") {
    params = corr_free;
    params_mean = params;
  }

  // Initializing theta: latent abilities
  arma::mat theta_prior_Sigma_chol = Cov_chol;
  theta_prior_Sigma_chol_inv = arma::inv(trimatl(theta_prior_Sigma_chol));
  theta_prior_Sigma_inv = theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;
  theta = arma::vec(n*m).fill(NA_REAL);

}

void Ifa::update_theta()
{
  if (model_type == "eifa" || model_type == "cifa" || model_type == "cifa_pred") {
    arma::mat aux_Sigma_inv_chol = arma::chol(LA.t()*LA + theta_prior_Sigma_inv, "lower");
    arma::mat aux_Sigma_chol = arma::inv(trimatl(aux_Sigma_inv_chol)).t();
    arma::mat aux_Sigma = aux_Sigma_chol * aux_Sigma_chol.t();
    arma::mat residual = vec2mat(z - arma::kron(c, ones_n), n, q);
    arma::vec theta_mu = arma::vectorise(residual * LA * aux_Sigma);
    theta = theta_mu + arma::kron(aux_Sigma_chol, eye_n) * arma::randn(n*m);
  }
}

void Ifa::update_c(const arma::vec& c_prior_mean, const arma::vec& c_prior_sd)
{
  arma::mat Theta = vec2mat(theta, n, m);
  arma::mat c_X = arma::kron(eye_q, ones_n);
  arma::vec c_Sigma_diag = 1 / (pow(c_prior_sd, -2) + n);
  arma::vec c_mu = c_X.t() * (z - arma::vectorise(Theta * LA.t()));
  c_mu %= c_Sigma_diag;
  c = c_mu + arma::randn(q) % sqrt(c_Sigma_diag);
}

void Ifa::update_a(const arma::vec& a_prior_mean, const arma::mat& A_prior_sd) {
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

void Ifa::update_cov_params(const double C, const double alpha, const double target,
    const double R_prior_eta, int index) {

  double accept;

  if (model_type == "cifa") {
    // Define current Sigma_proposal and propose new covariance matrix
    arma::mat Sigma_proposal = exp(logscale) * params_cov;
    arma::mat Sigma_proposal_chol = arma::chol(Sigma_proposal, "lower");
    arma::vec params_aux = params + Sigma_proposal_chol * arma::randn(params.size());
    arma::vec corr_free_aux = params_aux;
    arma::mat Corr_chol_aux = vec2chol_corr(params_aux, m);
    // Update covariance structure of latent factors (theta)
    arma::mat theta_prior_Sigma_chol_aux = Corr_chol_aux;
    theta_prior_Sigma_chol_aux.each_col() %= V_sd;
    arma::mat theta_prior_Sigma_chol_inv_aux = arma::inv(trimatl(theta_prior_Sigma_chol_aux));
    // Compute probability of acceptance for the RW update
    arma::mat Theta = vec2mat(theta, n, m);
    arma::mat Theta_mean = arma::zeros(n, m);
    accept = dmvnorm_cholinv(Theta.t(), Theta_mean.t(), theta_prior_Sigma_chol_inv_aux);
    accept += dlkj_corr_free2(corr_free_aux, m, R_prior_eta);
    accept -= dmvnorm_cholinv(Theta.t(), Theta_mean.t(), theta_prior_Sigma_chol_inv);
    accept -= dlkj_corr_free2(corr_free, m, R_prior_eta);
    // Accept or reject new proposal and update associated parameters
    if (accept > log(R::runif(0,1))) {
      params = params_aux;
      corr_free = corr_free_aux;
      Corr_chol = Corr_chol_aux;
      theta_prior_Sigma_chol_inv = theta_prior_Sigma_chol_inv_aux;
      theta_prior_Sigma_inv = theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;
    }
  }

  if (model_type != "eifa") {
    // Update scaling parameter
    double gamma = C / pow(index + 1, alpha);
    logscale += gamma * (std::min(exp(accept), 1.0) - target);
    // Update covariance of parameters with stochastic approximation
    arma::vec params_center = params - params_mean;
    params_mean += gamma * params_center;
    params_cov += gamma * (params_center * params_center.t() - params_cov);
  }

}

Rcpp::List Ifa::sample(
    arma::vec c_prior_mean, arma::vec c_prior_sd,
    arma::mat A_prior_mean, arma::mat A_prior_sd,
    int niter, int thin,
    double C, double alpha, double target, double R_prior_eta
    ) {

  int nsave = niter / thin;

  // Transformation of prior parameters
  arma::vec a_prior_mean = arma::vectorise(A_prior_mean.t());

  // Define matrices to save samples
  arma::mat c_samples(q, nsave);
  arma::mat a_samples(q*m, nsave);
  arma::mat theta_samples(n*m, nsave);
  arma::mat z_samples(q*n, nsave);
  arma::mat corr_samples(n_corr, nsave);


  for (int i = 0; i < niter; ++i) {
    // Update parameters
    update_theta();
    update_c(c_prior_mean, c_prior_sd);
    update_a(a_prior_mean, A_prior_sd);
    update_z();
    update_cov_params(C, alpha, target, R_prior_eta, i);
    // Save samples
    if (niter % thin == 0) {
      theta_samples.col(i/thin) = theta;
      c_samples.col(i/thin) = c;
      a_samples.col(i/thin) = arma::vectorise(LA);
      z_samples.col(i/thin) = z;
      corr_samples.col(i/thin) = trimatl2vec(Corr_chol * Corr_chol.t(), false);
    }
  }

  Rcpp::NumericMatrix theta_samples_rcpp = Rcpp::wrap(theta_samples.t());
  Rcpp::colnames(theta_samples_rcpp) = name_samples_mat(n, m, "Theta");
  Rcpp::NumericMatrix c_samples_rcpp = Rcpp::wrap(c_samples.t());
  Rcpp::colnames(c_samples_rcpp) = name_samples_vec(q, "c");
  Rcpp::NumericMatrix a_samples_rcpp = Rcpp::wrap(a_samples.t());
  Rcpp::colnames(a_samples_rcpp) = name_samples_mat(q, m, "A");
  Rcpp::NumericMatrix z_samples_rcpp = Rcpp::wrap(z_samples.t());
  Rcpp::colnames(z_samples_rcpp) = name_samples_mat(n, q, "Z");
  Rcpp::NumericMatrix corr_samples_rcpp = Rcpp::wrap(corr_samples.t());
  Rcpp::colnames(corr_samples_rcpp) = name_samples_lower(m, m, "Corr", false);

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("c") = c_samples_rcpp,
      Rcpp::Named("a") = a_samples_rcpp,
      Rcpp::Named("theta") = theta_samples_rcpp,
      Rcpp::Named("z") = z_samples_rcpp,
      Rcpp::Named("corr") = corr_samples_rcpp
      );

  Rcpp::StringVector myclass(2);
  myclass(0) = "spmirt.list";
  myclass(1) = "list";
  output.attr("class") = myclass;

  return output;
}

