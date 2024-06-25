
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
#include "name-samples.h"
#include "correlation.h"
#include "ifa.h"
#include "pdf.h"
// [[Rcpp::depends(RcppArmadillo)]]

Ifa::Ifa (Rcpp::NumericVector response, arma::mat predictors, arma::mat distances,
    int nobs, int nitems, int nfactors, int ngps,
    arma::mat constrain_L, arma::mat constrain_T, arma::vec constrain_V_sd,
    arma::mat adap_Sigma, double adap_scale,
    arma::vec c_ini, arma::mat A_ini, arma::mat R_ini,
    arma::mat B_ini, arma::vec sigmas_gp_ini, arma::vec phi_gp_ini,
    std::string mod_type):
  model_type(mod_type),
  y(response), dist(distances), X(predictors),
  n(nobs), q(nitems), m(nfactors), ngp(ngps), p(predictors.n_cols), n_corr((m-1)*m / 2),
  ones_n(arma::ones(n)),
  zeros_nm(arma::zeros(n*m)),
  eye_q(arma::eye(q,q)),
  eye_n(arma::eye(n,n)),
  eye_m(arma::eye(m,m)),
  low_thresh(Rcpp::NumericVector::create(R_NegInf, 0)),
  high_thresh(Rcpp::NumericVector::create(0, R_PosInf)),
  L(constrain_L),
  T(constrain_T),
  T_index(find(T != 0)),
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
    if (Rcpp::NumericVector::is_na(*ity)) {
      *it = arma::randn();
    } else {
      *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    }
    ity++;
  }
  // Initializing B: multivariate fixed effects
  B = B_ini;
  // Initializing correlation parameters of multivariate independent noise (V)
  // Rcpp::Rcout << "R_ini: "<< R_ini << std::endl;
  Corr_chol = arma::chol(R_ini, "lower");
  corr_free = chol_corr2vec(Corr_chol);
  // Rcpp::Rcout << "corr_free: "<< corr_free << std::endl;
  arma::mat Cov_chol = Corr_chol;
  Cov_chol.each_col() %= V_sd;
  // Initializing sigmas and phis: spatial covariance parameters of MGP
  mgp_sd = sigmas_gp_ini;
  mgp_phi = phi_gp_ini;
  // Update full covariance structure of latent abilities depending of model type
  if (model_type == "spifa" || model_type == "spifa_pred") {
    arma::mat mgp_Sigma(n*ngp, n*ngp, arma::fill::zeros);
    for (int i = 0; i < ngp; ++i) {
      mgp_Sigma.submat(i*n, i*n, arma::size(n, n)) = exp(- dist / mgp_phi(i));
    }
    T(T_index) = mgp_sd;
    mgp_Sigma = TST(mgp_Sigma, T);
    arma::mat theta_prior_Sigma = mgp_Sigma + arma::kron(Cov_chol * Cov_chol.t(), eye_n);
    arma::mat theta_prior_Sigma_chol = arma::chol(theta_prior_Sigma, "lower");
    theta_prior_Sigma_chol_inv = arma::inv(trimatl(theta_prior_Sigma_chol));
    theta_prior_Sigma_inv =
      theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;
    // Rcpp::Rcout << theta_prior_Sigma.diag() << std::endl;
  } else if (model_type == "eifa" || model_type == "cifa" || model_type == "cifa_pred") {

    arma::mat theta_prior_Sigma_chol = Cov_chol;
    theta_prior_Sigma_chol_inv = arma::inv(trimatl(theta_prior_Sigma_chol));
    theta_prior_Sigma_inv = theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;
    // Rcpp::Rcout << "chol: "<< Cov_chol << std::endl;
    // Rcpp::Rcout << "sigmachol: "<< theta_prior_Sigma_chol << std::endl;
    // Rcpp::Rcout << "sigmainv: "<< theta_prior_Sigma_inv << std::endl;

  }
  // Initializing mean of latent abilities depending of model type
  if (model_type == "cifa_pred" || model_type == "spifa_pred") {
    theta_prior_mean = arma::vectorise(X * B);
  } else {
    theta_prior_mean = arma::zeros(n*m);
  }
  // Initializing theta: latent abilities
  theta = arma::vec(n*m).fill(NA_REAL);

  // Join covariance parameters depending of model type
  if (model_type == "spifa" || model_type == "spifa_pred") {
    params = arma::join_cols(log(mgp_sd), log(mgp_phi));
    params = arma::join_cols(params, corr_free);
    // Rcpp::Rcout << params << std::endl;
  } else if (model_type == "eifa" || model_type == "cifa" || model_type == "cifa_pred") {
    params = corr_free;
  }
  // Variables for adaptive MH
  params_mean = params;
  // Rcpp::Rcout << "corr_free: "<< params << std::endl;
}

Ifa::Ifa(Rcpp::NumericVector response, arma::mat predictors, arma::mat distances,
      int nobs, int nitems, int nfactors, int ngps,
      arma::mat constrain_L, arma::mat constrain_T, arma::mat constrain_V_sd,
      std::string mod_type):
  model_type(mod_type),
  y(response), dist(distances), X(predictors),
  n(nobs), q(nitems), m(nfactors), ngp(ngps), p(predictors.n_cols), n_corr((m-1)*m / 2),
  ones_n(arma::ones(n)),
  zeros_nm(arma::zeros(n*m)),
  eye_q(arma::eye(q,q)),
  eye_n(arma::eye(n,n)),
  eye_m(arma::eye(m,m)),
  low_thresh(Rcpp::NumericVector::create(R_NegInf, 0)),
  high_thresh(Rcpp::NumericVector::create(0, R_PosInf)),
  L(constrain_L),
  T(constrain_T),
  T_index(find(T != 0)),
  V_sd(constrain_V_sd)
{
}



void Ifa::update_theta()
{
  if (model_type == "eifa" || model_type == "cifa" || model_type == "cifa_pred") {
    arma::mat aux_Sigma_inv_chol = arma::chol(LA.t()*LA + theta_prior_Sigma_inv, "lower");
    arma::mat aux_Sigma_chol = arma::inv(trimatl(aux_Sigma_inv_chol)).t();
    arma::mat aux_Sigma = aux_Sigma_chol * aux_Sigma_chol.t();
    arma::mat residual = vec2mat(z - arma::kron(c, ones_n), n, q);
    arma::vec theta_mu = arma::vectorise(residual * LA * aux_Sigma);
    if (model_type == "cifa_pred") {
      theta_mu += arma::vectorise(X * B * theta_prior_Sigma_inv * aux_Sigma);
    }
    theta = theta_mu + arma::kron(aux_Sigma_chol, eye_n) * arma::randn(n*m);
  } else if (model_type == "spifa" || model_type == "spifa_pred") {
    arma::mat theta_Sigma_inv = arma::kron(LA.t() * LA, eye_n) +  theta_prior_Sigma_inv;
    arma::mat theta_Sigma_inv_chol = arma::chol(theta_Sigma_inv, "lower");
    arma::mat theta_Sigma_chol = arma::inv(trimatl(theta_Sigma_inv_chol)).t();
    arma::vec theta_mu = arma::kron(LA.t(), eye_n) * (z - arma::kron(c, ones_n));
    if (model_type == "spifa_pred") {
      theta_mu += theta_prior_Sigma_inv * theta_prior_mean;
    }
    theta_mu = theta_Sigma_chol * theta_Sigma_chol.t() * theta_mu;
    theta = theta_mu +  theta_Sigma_chol * arma::randn(n*m);
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

void Ifa::update_z()
{
  arma::mat Theta = vec2mat(theta, n, m);
  arma::vec z_mu = arma::kron(c, ones_n) + arma::vectorise(Theta * LA.t());
  Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itz_mu = z_mu.begin();
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    if (Rcpp::NumericVector::is_na(*ity)) {
      *it = *itz_mu + arma::randn();
    } else {
      *it = RcppTN::rtn1(*itz_mu, 1.0, low_thresh[*ity], high_thresh[*ity]);
    }
    ity++;
    itz_mu++;
  }
}

void Ifa::update_B(const arma::mat& B_prior_sd)
{
    arma::mat Theta = vec2mat(theta, n, m);
  if (model_type == "cifa_pred") {
    arma::mat B_Sigma_inv = arma::kron(theta_prior_Sigma_inv, X.t()*X);
    B_Sigma_inv.diag() += arma::vectorise(arma::square(B_prior_sd));
    arma::mat B_Sigma_inv_chol = arma::chol(B_Sigma_inv, "lower");
    arma::mat B_Sigma_chol = arma::inv(trimatl(B_Sigma_inv_chol)).t();
    arma::mat B_Sigma = B_Sigma_chol * B_Sigma_chol.t();
    arma::vec B_mean = B_Sigma * arma::vectorise(X.t() * Theta * theta_prior_Sigma_inv);
    arma::vec betas = B_mean + B_Sigma_chol * arma::randn(p*m);
    B = vec2mat(betas, p, m);
  } else if (model_type == "spifa_pred") {
    arma::mat Xt_kron = arma::kron(eye_m, X.t());
    arma::mat B_Sigma_inv = Xt_kron * theta_prior_Sigma_inv * Xt_kron.t();
    B_Sigma_inv.diag() += arma::vectorise(arma::square(B_prior_sd));
    arma::mat B_Sigma_inv_chol = arma::chol(B_Sigma_inv, "lower");
    arma::mat B_Sigma_chol = arma::inv(trimatl(B_Sigma_inv_chol)).t();
    arma::mat B_Sigma = B_Sigma_chol * B_Sigma_chol.t();
    arma::vec B_mean = arma::vectorise(X.t()*vec2mat(theta_prior_Sigma_inv*theta, n, m));
    B_mean = B_Sigma * B_mean;
    arma::vec betas = B_mean + B_Sigma_chol * arma::randn(p*m);
    B = vec2mat(betas, p, m);
  }
}

void Ifa::update_cov_params(
    const arma::vec sigmas_gp_mean, const arma::vec sigmas_gp_sd,
    const arma::vec phi_gp_mean, const arma::vec phi_gp_sd,
    const double R_prior_eta,
    const double C, const double alpha, const double target, int index)
{
  double accept;

  if (model_type == "cifa" || model_type == "cifa_pred") {
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
    arma::mat Theta_mean = vec2mat(theta_prior_mean, n, m);
    accept = dmvnorm_cholinv(Theta.t(), Theta_mean.t(), theta_prior_Sigma_chol_inv_aux);
    accept += dlkj_corr_free2(corr_free_aux, m, R_prior_eta);
    accept -= dmvnorm_cholinv(Theta.t(), Theta_mean.t(), theta_prior_Sigma_chol_inv);
    accept -= dlkj_corr_free2(corr_free, m, R_prior_eta);
    // Accept or reject new proposal and update associated parameters
    if (accept > log(R::runif(0,1))) {
      // All parameters
      params = params_aux;
      // Correlation updates
      corr_free = corr_free_aux;
      Corr_chol = Corr_chol_aux;
      // Latent abilities prior update
      theta_prior_Sigma_chol_inv = theta_prior_Sigma_chol_inv_aux;
      theta_prior_Sigma_inv = theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;
    }
  } else if (model_type == "spifa" || model_type == "spifa_pred") {
    const int nsig = mgp_sd.n_elem;
    // Define current Sigma_proposal and update covariance matrix
    arma::mat Sigma_proposal = exp(logscale) * params_cov;
    arma::mat Sigma_proposal_chol = arma::chol(Sigma_proposal, "lower");
    arma::vec params_aux = params + Sigma_proposal_chol * arma::randn(params.size());
    // Build proposed mgp covariance structure
    arma::vec mgp_sd_aux = exp(params_aux.subvec(0, arma::size(nsig,1)));
    arma::vec mgp_phi_aux = exp(params_aux.subvec(nsig, arma::size(ngp,1)));
    arma::mat mgp_Sigma_aux(n*ngp, n*ngp, arma::fill::zeros);
    for (int j = 0; j < ngp; ++j) {
      mgp_Sigma_aux.submat(j*n, j*n, arma::size(n, n)) = exp(-dist / mgp_phi_aux(j));
    }
    arma::mat T_aux = arma::zeros(m, ngp);
    T_aux(T_index) = mgp_sd_aux;
    mgp_Sigma_aux = TST(mgp_Sigma_aux, T_aux);
    // Build proposed multivariate noise covariance structure
    arma::vec corr_free_aux = params_aux.subvec(nsig + ngp, params.size()-1);
    arma::mat Corr_chol_aux = vec2chol_corr(corr_free_aux, m);
    arma::mat Cov_chol_aux = Corr_chol_aux;
    Cov_chol_aux.each_col() %= V_sd;
    arma::mat Cov_aux = Cov_chol_aux * Cov_chol_aux.t();
    // Build full proposed covariance structure
    arma::mat theta_prior_Sigma_aux = mgp_Sigma_aux + arma::kron(Cov_aux, eye_n);
    arma::mat theta_prior_Sigma_chol_aux = arma::chol(theta_prior_Sigma_aux, "lower");
    arma::mat theta_prior_Sigma_chol_inv_aux =
      arma::inv(trimatl(theta_prior_Sigma_chol_aux));
    // Update
    accept = dmvnorm_cholinv(theta, theta_prior_mean, theta_prior_Sigma_chol_inv_aux, true);
    accept += dlkj_corr_free2(corr_free_aux, m, R_prior_eta, true);
    accept -= dmvnorm_cholinv(theta, theta_prior_mean, theta_prior_Sigma_chol_inv, true);
    accept -= dlkj_corr_free2(corr_free, m, R_prior_eta, true);
    for (int k = 0; k < nsig; ++k) {
     accept += R::dnorm(log(mgp_sd_aux(k)), log(sigmas_gp_mean(k)), sigmas_gp_sd(k), true);
     accept -= R::dnorm(log(mgp_sd(k)), log(sigmas_gp_mean(k)), sigmas_gp_sd(k), true);
    }
    for (int j = 0; j < ngp; ++j) {
     accept += R::dnorm(log(mgp_phi_aux(j)), log(phi_gp_mean(j)), phi_gp_sd(j), true);
     accept -= R::dnorm(log(mgp_phi(j)), log(phi_gp_mean(j)), phi_gp_sd(j), true);
    }
    // Rcpp::Rcout << accept << std::endl;
    if (accept > log(R::runif(0,1))) {
      // All parameters
      params = params_aux;
      // Correlation updates
      corr_free = corr_free_aux;
      Corr_chol = Corr_chol_aux;
      // Gaussian processes parameters updates
      mgp_sd = mgp_sd_aux;
      T = T_aux;
      mgp_phi = mgp_phi_aux;
      // Latent abilities prior update
      theta_prior_Sigma_chol_inv = theta_prior_Sigma_chol_inv_aux;
      theta_prior_Sigma_inv =
        theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;
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
    int niter, int thin, bool standardize,
    arma::vec c_prior_mean, arma::vec c_prior_sd,
    arma::mat A_prior_mean, arma::mat A_prior_sd,
    double R_prior_eta,
    arma::mat B_prior_mean, arma::mat B_prior_sd,
    arma::vec sigmas_gp_mean, arma::vec sigmas_gp_sd,
    arma::vec phi_gp_mean, arma::vec phi_gp_sd,
    double C, double alpha, double target)
{

  int nsave = niter / thin;

  // Transformation of prior parameters
  arma::vec a_prior_mean = arma::vectorise(A_prior_mean.t());

  int n_chol_corr = m*(m+1)/2;
  // Define matrices to save samples
  arma::mat c_samples(q, nsave);
  arma::mat a_samples(q*m, nsave);
  arma::mat theta_samples(n*m, nsave);
  arma::mat z_samples(q*n, nsave);
  arma::mat corr_samples(n_corr, nsave);
  arma::mat corr_chol_samples(n_chol_corr, nsave);
  arma::mat mgp_sd_samples(mgp_sd.n_elem, nsave);
  arma::mat mgp_phi_samples(ngp, nsave);
  arma::mat betas_samples(p*m, nsave);

  for (int i = 0; i < niter; ++i) {
    // Update parameters
    update_theta();
    update_c(c_prior_mean, c_prior_sd);
    update_a(a_prior_mean, A_prior_sd);
    update_z();
    update_B(B_prior_sd);
    update_cov_params(sigmas_gp_mean, sigmas_gp_sd, phi_gp_mean, phi_gp_sd,
        R_prior_eta, C, alpha, target, i);

    // Save samples
    if (niter % thin == 0) {
      theta_samples.col(i/thin) = theta;
      c_samples.col(i/thin) = c;
      a_samples.col(i/thin) = arma::vectorise(LA);
      z_samples.col(i/thin) = z;
      corr_chol_samples.col(i/thin) = trimatl2vec(Corr_chol, true);
      corr_samples.col(i/thin) = trimatl2vec(Corr_chol * Corr_chol.t(), false);
      mgp_sd_samples.col(i/thin) = mgp_sd;
      mgp_phi_samples.col(i/thin) = mgp_phi;
      betas_samples.col(i/thin) = arma::vectorise(B);
    }
  }

  if (standardize && model_type != "eifa" && model_type != "cifa") {
    arma::vec theta_samples_sd(m, arma::fill::zeros);
    arma::umat T_sub = arma::ind2sub(size(T), T_index); // indices to row-column index
    Rcpp::Rcout << "Standardixing" << std::endl;
    for (int i = 0; i < m; ++i) {
      // get variance of the latent abilities
      arma::vec sds =
        arma::stddev(theta_samples.submat(i*n, nsave/2, arma::size(n, nsave/2))).t();
      theta_samples_sd(i) = arma::mean(sds);
      //  standardize latent abilities
      theta_samples.submat(i*n, 0, arma::size(n, nsave)) /= theta_samples_sd(i);
      // standardize discrimation parameters
      a_samples.submat(i*q, 0, arma::size(q, nsave)) *= theta_samples_sd(i);
      // standardize betas
      betas_samples.submat(i*p, 0, arma::size(p, nsave)) /= theta_samples_sd(i);
      // standardize gaussian process variance parameter if diagonal
      mgp_sd_samples.rows(find(T_sub.row(0) == i)) /= theta_samples_sd(i);
    }
    // new standard deviation for the multivariate residual term
    V_sd = V_sd / theta_samples_sd;
  }

  Rcpp::NumericMatrix theta_samples_rcpp = Rcpp::wrap(theta_samples.t());
  Rcpp::colnames(theta_samples_rcpp) = name_samples_mat(n, m, "Theta");
  Rcpp::NumericMatrix c_samples_rcpp = Rcpp::wrap(c_samples.t());
  Rcpp::colnames(c_samples_rcpp) = name_samples_vec(q, "c");
  Rcpp::NumericMatrix a_samples_rcpp = Rcpp::wrap(a_samples.t());
  Rcpp::colnames(a_samples_rcpp) = name_samples_mat(q, m, "A");
  Rcpp::NumericMatrix z_samples_rcpp = Rcpp::wrap(z_samples.t());
  Rcpp::colnames(z_samples_rcpp) = name_samples_mat(n, q, "Z");
  Rcpp::NumericMatrix corr_chol_samples_rcpp = Rcpp::wrap(corr_chol_samples.t());
  Rcpp::colnames(corr_chol_samples_rcpp) = name_samples_lower(m, m, "Corr_chol", true);
  Rcpp::NumericMatrix corr_samples_rcpp = Rcpp::wrap(corr_samples.t());
  Rcpp::colnames(corr_samples_rcpp) = name_samples_lower(m, m, "Corr", false);
  Rcpp::NumericMatrix mgp_sd_samples_rcpp = Rcpp::wrap(mgp_sd_samples.t());
  Rcpp::colnames(mgp_sd_samples_rcpp) =
    name_samples_mat(m, ngp, "T")[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(T_index))];
  Rcpp::NumericMatrix mgp_phi_samples_rcpp = Rcpp::wrap(mgp_phi_samples.t());
  Rcpp::colnames(mgp_phi_samples_rcpp) = name_samples_vec(ngp, "mgp_phi");
  Rcpp::NumericMatrix betas_samples_rcpp = Rcpp::wrap(betas_samples.t());
  Rcpp::colnames(betas_samples_rcpp) = name_samples_mat(p, m, "B");

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("c") = c_samples_rcpp,
      Rcpp::Named("a") = a_samples_rcpp,
      Rcpp::Named("theta") = theta_samples_rcpp,
      Rcpp::Named("z") = z_samples_rcpp,
      Rcpp::Named("corr_chol") = corr_chol_samples_rcpp,
      Rcpp::Named("corr") = corr_samples_rcpp,
      Rcpp::Named("mgp_sd") = mgp_sd_samples_rcpp,
      Rcpp::Named("mgp_phi") = mgp_phi_samples_rcpp,
      Rcpp::Named("betas") = betas_samples_rcpp
      );

  Rcpp::StringVector myclass(2);
  myclass(0) = "spifa.list";
  myclass(1) = "list";
  output.attr("class") = myclass;
  output.attr("V_sd") = V_sd;

  return output;
}

Rcpp::List Ifa::predict(arma::mat samples_theta, arma::mat samples_corr_chol,
    arma::mat samples_corr,
    arma::mat samples_mgp_sd, arma::mat samples_mgp_phi, arma::mat samples_betas,
    arma::mat newpredictors, arma::mat newdist, arma::mat cross_distances,
    int npred, int niter, int burnin, int thin)
{

  // dimensions
  int nsave = (niter - burnin) / thin;
  arma::mat prediction(npred * m, nsave);

  int isave = 0;

  if (model_type == "cifa_pred") {
    for (int i = burnin; i < niter; i += thin) {
      // compute the mean
      arma::mat B = vec2mat(samples_betas.col(i), p, m);
      arma::vec pred_mean = arma::vectorise(newpredictors * B);
      // build auxiliary cholesky of covariance
      arma::mat Corr_chol_aux = vec2trimatl(samples_corr_chol.col(i), m);
      arma::mat pred_cov_chol_aux = Corr_chol_aux;
      pred_cov_chol_aux.each_col() %= V_sd;
      // samples prediction
      prediction.col(isave) =
        pred_mean + arma::vectorise(arma::randn(npred, m) * pred_cov_chol_aux.t());
      isave++;
    }
  }

  if (model_type == "spifa" || model_type == "spifa_pred") {
    arma::mat eye_npred = arma::eye(npred, npred);
    for (int i = burnin; i < niter; i += thin) {
      // compute means
      arma::vec obs_mean(n*m, arma::fill::zeros);
      arma::vec new_mean(npred*m, arma::fill::zeros);
      if (model_type == "spifa_pred") {
        arma::mat B = vec2mat(samples_betas.col(i), p, m);
        obs_mean = arma::vectorise(X * B);
        new_mean = arma::vectorise(newpredictors * B);
      }
      // compute correlation of residual term
      arma::mat Cov_aux = vec2trimatl(samples_corr.col(i), m, false);
      Cov_aux.diag().ones();
      Cov_aux.each_col() %= V_sd;
      Cov_aux.each_row() %= V_sd.t();
      // compute variances and covariances
      arma::mat obs_mgp_Sigma(n*ngp, n*ngp, arma::fill::zeros);
      arma::mat new_mgp_Sigma(npred*ngp, npred*ngp, arma::fill::zeros);
      arma::mat new_mgp_Cov(npred*ngp, n*ngp, arma::fill::zeros);
      for (int j = 0; j < ngp; ++j) {
        obs_mgp_Sigma.submat(j*n, j*n, arma::size(n, n)) =
          exp(- dist / samples_mgp_phi(j,i));
        new_mgp_Sigma.submat(j*npred, j*npred, arma::size(npred, npred)) =
          exp(- newdist / samples_mgp_phi(j,i));
        new_mgp_Cov.submat(j*npred, j*n, arma::size(npred, n)) =
          exp(- cross_distances / samples_mgp_phi(j,i));
      }
      T(T_index) = samples_mgp_sd.col(i);
      obs_mgp_Sigma = TST(obs_mgp_Sigma, T) + arma::kron(Cov_aux, eye_n);
      new_mgp_Sigma = TST(new_mgp_Sigma, T) + arma::kron(Cov_aux, eye_npred);
      new_mgp_Cov = arma::kron(T, eye_npred) * new_mgp_Cov * arma::kron(T.t(), eye_n);
      // inv(chol(S_22)) * S_21
      arma::mat obs_mgp_Sigma_chol_inv =
        arma::inv(trimatl(arma::chol(obs_mgp_Sigma, "lower")));
      arma::mat L_inv_S21 = obs_mgp_Sigma_chol_inv * new_mgp_Cov.t();
      arma::mat L_inv_res = obs_mgp_Sigma_chol_inv * (samples_theta.col(i) - obs_mean);
      // compute prediction mean and variance
      arma::vec pred_mean = new_mean + L_inv_S21.t() * L_inv_res;
      arma::mat pred_var = new_mgp_Sigma - L_inv_S21.t() * L_inv_S21;
      // sample prediction
      if (true) {
        prediction.col(isave) =
          pred_mean + sqrt(pred_var.diag()) % arma::randn(npred * m);
      } else {
        prediction.col(isave) =
          pred_mean + arma::chol(pred_var, "lower") * arma::randn(npred * m);
      }
      isave++;
    }
  }

  Rcpp::NumericMatrix prediction_rcpp = Rcpp::wrap(prediction.t());
  Rcpp::colnames(prediction_rcpp) = name_samples_mat(npred, m, "Theta");

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("theta") = prediction_rcpp
      );

  return output;
}

Rcpp::List Ifa::predict2(arma::mat samples_theta, arma::mat samples_corr_chol,
    arma::mat samples_corr,
    arma::mat samples_mgp_sd, arma::mat samples_mgp_phi, arma::mat samples_betas,
    arma::mat newpredictors, arma::mat newdist, arma::mat cross_distances,
    int npred, int niter, int burnin, int thin)
{

  // dimensions
  int nsave = (niter - burnin) / thin;
  arma::mat prediction(npred * m, nsave);

  int isave = 0;

  if (model_type == "cifa_pred") {
    for (int i = burnin; i < niter; i += thin) {
      // compute the mean
      arma::mat B = vec2mat(samples_betas.col(i), p, m);
      arma::vec pred_mean = arma::vectorise(newpredictors * B);
      // build auxiliary cholesky of covariance
      arma::mat Corr_chol_aux = vec2trimatl(samples_corr_chol.col(i), m);
      arma::mat pred_cov_chol_aux = Corr_chol_aux;
      pred_cov_chol_aux.each_col() %= V_sd;
      // samples prediction
      prediction.col(isave) =
        pred_mean + arma::vectorise(arma::randn(npred, m) * pred_cov_chol_aux.t());
      isave++;
    }
  }

  if (model_type == "spifa" || model_type == "spifa_pred") {
    arma::mat eye_npred = arma::eye(npred, npred);
    for (int i = burnin; i < niter; i += thin) {
      // compute means
      arma::vec obs_mean(n*m, arma::fill::zeros);
      arma::vec new_mean(npred*m, arma::fill::zeros);
      if (model_type == "spifa_pred") {
        arma::mat B = vec2mat(samples_betas.col(i), p, m);
        obs_mean = arma::vectorise(X * B);
        // new_mean = arma::vectorise(newpredictors * B);
      }
      // compute correlation of residual term
      arma::mat Cov_aux = vec2trimatl(samples_corr.col(i), m, false);
      Cov_aux.diag().ones();
      Cov_aux.each_col() %= V_sd;
      Cov_aux.each_row() %= V_sd.t();
      // compute variances and covariances
      arma::mat obs_mgp_Sigma(n*ngp, n*ngp, arma::fill::zeros);
      arma::mat new_mgp_Sigma(npred*ngp, npred*ngp, arma::fill::zeros);
      arma::mat new_mgp_Cov(npred*ngp, n*ngp, arma::fill::zeros);
      for (int j = 0; j < ngp; ++j) {
        obs_mgp_Sigma.submat(j*n, j*n, arma::size(n, n)) =
          exp(- dist / samples_mgp_phi(j,i));
        new_mgp_Sigma.submat(j*npred, j*npred, arma::size(npred, npred)) =
          exp(- newdist / samples_mgp_phi(j,i));
        new_mgp_Cov.submat(j*npred, j*n, arma::size(npred, n)) =
          exp(- cross_distances / samples_mgp_phi(j,i));
      }
      T(T_index) = samples_mgp_sd.col(i);
      obs_mgp_Sigma = TST(obs_mgp_Sigma, T) + arma::kron(Cov_aux, eye_n);
      new_mgp_Sigma = TST(new_mgp_Sigma, T);
      new_mgp_Cov = arma::kron(T, eye_npred) * new_mgp_Cov * arma::kron(T.t(), eye_n);
      // inv(chol(S_22)) * S_21
      arma::mat obs_mgp_Sigma_chol_inv =
        arma::inv(trimatl(arma::chol(obs_mgp_Sigma, "lower")));
      arma::mat L_inv_S21 = obs_mgp_Sigma_chol_inv * new_mgp_Cov.t();
      arma::mat L_inv_res = obs_mgp_Sigma_chol_inv * (samples_theta.col(i) - obs_mean);
      // compute prediction mean and variance
      arma::vec pred_mean = new_mean + L_inv_S21.t() * L_inv_res;
      arma::mat pred_var = new_mgp_Sigma - L_inv_S21.t() * L_inv_S21;
      // sample prediction
      if (true) {
        prediction.col(isave) =
          pred_mean + sqrt(pred_var.diag()) % arma::randn(npred * m);
      } else {
        prediction.col(isave) =
          pred_mean + arma::chol(pred_var, "lower") * arma::randn(npred * m);
      }
      isave++;
    }
  }

  Rcpp::NumericMatrix prediction_rcpp = Rcpp::wrap(prediction.t());
  Rcpp::colnames(prediction_rcpp) = name_samples_mat(npred, m, "Theta");

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("theta") = prediction_rcpp
      );

  return output;
}

