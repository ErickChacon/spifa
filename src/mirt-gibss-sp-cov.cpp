
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
#include "pdf.h"
#include "correlation.h"
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Spatial Confirmatory Item Factor Analysis
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
Rcpp::List ifa_gibbs_sp_cov(Rcpp::NumericVector y, arma::mat X, arma::mat dist,
    int n, int q, int m,
    arma::vec mgp_phi, arma::vec mgp_sd, arma::mat Corr, double sd_fix,
    arma::mat sigma_prop, arma::mat L, arma::mat T,
    double target = 0.234, int niter = 1000) {

  // Constanst objects
  arma::vec ones_n(n, arma::fill::ones);
  arma::mat eye_q = arma::eye<arma::mat>(q,q);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  arma::mat eye_m = arma::eye<arma::mat>(m,m);
  arma::vec zeros_nm(n*m, arma::fill::zeros);
  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);
  Rcpp::NumericVector lows = low_thresh[y];
  Rcpp::NumericVector highs = high_thresh[y];
  const int n_corr = round((m-1) * m / 2);
  const int p = X.n_cols;

  // Hyperparameters of priors
  double prior_corr_eta = 1.5;

  // Initializing c: difficulty parameters
  arma::vec c = arma::randn(q) * 0.3;
  arma::vec c_mu(q);
  arma::mat c_X = arma::kron(eye_q, ones_n);
  arma::mat c_samples(q, niter);

  // Initializing a: discrimation parameters
  // arma::mat A = arma::randn(q, m) * 0.5;
  arma::mat A = arma::ones(q, m) * 0.5;
  arma::mat LA = L % A;
  arma::vec a = arma::vectorise(A.t());
  arma::vec a_mu(q*m);
  arma::mat a_Sigma(q*m, q*m, arma::fill::zeros);
  arma::mat a_Sigma_chol(q*m, q*m, arma::fill::zeros);
  arma::mat a_samples(q*m, niter);

  // Initializing z: auxiliary variables
  arma::vec z = arma::zeros<arma::vec>(y.size());
  arma::vec z_mu(y.size());
  Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itz_mu;
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  arma::mat z_samples(q*n, niter);

  // Initialize beta: fixed effect parameters
  arma::vec beta = arma::randn(p * m) * 0.1;
  arma::mat Beta = vec2mat(beta, p, m);
  arma::mat beta_samples(p * m, niter);

  // Initializing sigmas and phis: spatial covariance parameters of GP
  arma::vec mgp_sd_log = log(mgp_sd);
  arma::vec mgp_phi_log = log(mgp_phi);
  arma::mat mgp_Sigma(n*m, n*m, arma::fill::zeros);
  for (int i = 0; i < m; ++i) {
    mgp_Sigma.submat(i*n, i*n, arma::size(n, n)) =
      pow(mgp_sd(i), 2) * exp(- dist / mgp_phi(i));
  }
  arma::mat mgp_sd_samples(m, niter);
  arma::mat mgp_phi_samples(m, niter);


  // Initializing correlation parameters: correlation of multivariate noise
  // arma::mat Cov = Corr * pow(sd_fix, 2);
  arma::vec corr_free(n_corr, arma::fill::randn);
  arma::mat Corr_chol = vec2chol_corr(corr_free, m);
  arma::mat Cov_chol = Corr_chol * sd_fix;
  arma::mat Cov = Cov_chol * Cov_chol.t();
  // arma::mat Cov_chol = Corr_chol;
  // Cov_chol.each_col() %= sigmas;
  arma::mat Cov_chol_inv = arma::inv(trimatl(Cov_chol));
  arma::mat Cov_inv = Cov_chol_inv.t() * Cov_chol_inv;
  arma::mat corr_samples(n_corr, niter);
  // arma::mat corr_samples2(n_corr, iter);
  // double corr_free_logscale = 0;


  // Update covariance structure
  arma::mat theta_prior_Sigma = mgp_Sigma + arma::kron(Cov, eye_n);
  arma::mat theta_prior_Sigma_chol = arma::chol(theta_prior_Sigma, "lower");
  arma::mat theta_prior_Sigma_chol_inv = arma::inv(trimatl(theta_prior_Sigma_chol));
  arma::mat theta_prior_Sigma_inv =
    theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;

  // Set parameters for theta: latent abilities
  arma::vec theta(n*m);
  arma::mat Theta(n, m);
  arma::mat theta_Sigma_inv(n*m, n*m);
  arma::mat theta_Sigma_chol(n*m, n*m);
  arma::vec theta_mu(n*m);
  arma::mat theta_samples(n*m, niter);


  // Initializing parameter variables for Metropolis Hastings
  arma::vec params = arma::join_cols(mgp_sd_log, mgp_phi_log);
  params = arma::join_cols(params, corr_free);
  arma::vec params_mean = params;
  arma::mat params_cov = sigma_prop;
  double logscale = 0;

  // Adaptive stepsize
  // As n increases alpha need to go towards 1
  // For n=35, alpha = 0.8 and C = 0.7
  const double alpha = 0.8;
  const double C = 0.7;
  double gamma;
  double accept;


  for (int i = 0; i < niter; ++i) {

    // Updating latent habilities (theta)
    theta_Sigma_inv = arma::kron(LA.t() * LA, eye_n) +  theta_prior_Sigma_inv;
    arma::mat theta_Sigma_inv_chol = arma::chol(theta_Sigma_inv, "lower");
    theta_Sigma_chol = arma::inv(trimatl(theta_Sigma_inv_chol)).t();
    theta_mu = arma::kron(LA.t(), eye_n) * (z - arma::kron(c, ones_n));
    theta_mu = theta_Sigma_chol * theta_Sigma_chol.t() * theta_mu;
    theta = theta_mu +  theta_Sigma_chol * arma::randn<arma::vec>(n*m);
    Theta = vec2mat(theta, n, m);

    // Updating difficulty parameters (c)
    c_mu = c_X.t() * (z - arma::vectorise(Theta * LA.t())) / (n + 1);
    c = c_mu + arma::randn<arma::vec>(q) / sqrt(n + 1);

    // Updating discrimation parameters (a)
    arma::mat a_Sigma_Lt(q*m, q*m, arma::fill::zeros);
    arma::mat cross_Theta = Theta.t() * Theta;
    for (int j = 0; j < q; ++j) {
      arma::mat aux_a_Sigma_inv = cross_Theta;
      aux_a_Sigma_inv.rows(find(L.row(j) == 0)) *= 0;
      aux_a_Sigma_inv.cols(find(L.row(j) == 0)) *= 0;
      aux_a_Sigma_inv.diag() += 1 / 0.5;
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
    a_mu = a_Sigma_Lt * arma::kron(eye_q, Theta.t()) * (z - arma::kron(c, ones_n));
    a = a_mu + a_Sigma_chol * arma::randn<arma::vec>(q*m);
    A = vec2matt(a, m, q);
    LA = L % A;

    // Auxiliary variables (z)
    z_mu = arma::kron(c, ones_n) + arma::vectorise(Theta * LA.t());
    ity = y.begin();
    itz_mu = z_mu.begin();
    for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
      *it = RcppTN::rtn1(*itz_mu, 1.0, low_thresh[*ity], high_thresh[*ity]);
      ity++;
      itz_mu++;
    }

    // Define current Sigma_proposal and update covariance matrix
    arma::mat Sigma_proposal = exp(logscale) * params_cov;
    arma::mat Sigma_proposal_chol = arma::chol(Sigma_proposal, "lower");
    arma::vec params_aux = params +
      Sigma_proposal_chol * arma::randn<arma::vec>(params.size());
    arma::vec mgp_sd_aux = exp(params_aux.subvec(0, m-1));
    arma::vec mgp_phi_aux = exp(params_aux.subvec(m, 2*m-1));
    arma::vec corr_free_aux = params_aux.subvec(2*m, params.size()-1);
    arma::mat Corr_chol_aux = vec2chol_corr(corr_free_aux, m);
    arma::mat Cov_chol_aux = Corr_chol_aux * sd_fix;
    arma::mat Cov_aux = Cov_chol_aux * Cov_chol_aux.t();
    arma::mat mgp_Sigma_aux(n*m, n*m, arma::fill::zeros);
    for (int j = 0; j < m; ++j) {
      mgp_Sigma_aux.submat(j*n, j*n, arma::size(n, n)) =
        pow(mgp_sd_aux(j), 2) * exp(- dist / mgp_phi_aux(j));
    }
    // Update covariance structure
    arma::mat theta_prior_Sigma_aux = mgp_Sigma_aux + arma::kron(Cov_aux, eye_n);
    arma::mat theta_prior_Sigma_chol_aux = arma::chol(theta_prior_Sigma_aux, "lower");
    arma::mat theta_prior_Sigma_chol_inv_aux =
      arma::inv(trimatl(theta_prior_Sigma_chol_aux));
    // Update
    accept = dmvnorm_cholinv(theta, zeros_nm, theta_prior_Sigma_chol_inv_aux, true);
    accept += dlkj_corr_free(corr_free_aux, m, prior_corr_eta, true);
    accept -= dmvnorm_cholinv(theta, zeros_nm, theta_prior_Sigma_chol_inv, true);
    accept -= dlkj_corr_free(corr_free, m, prior_corr_eta, true);
    for (int j = 0; j < m; ++j) {
     accept += R::dnorm(log(mgp_sd_aux(j)), log(pow(0.6,2)), 0.2, true);
     accept += R::dnorm(log(mgp_phi_aux(j)), log(0.04), 0.2, true);
     accept -= R::dnorm(log(mgp_sd(j)), log(pow(0.6,2)), 0.2, true);
     accept -= R::dnorm(log(mgp_phi(j)), log(0.04), 0.2, true);
    }
    if (accept > log(R::runif(0,1))) {
      // Correlation updates
      corr_free = corr_free_aux;
      Corr_chol = Corr_chol_aux;
      Cov_chol = Cov_chol_aux;
      Cov_chol_inv = arma::inv(trimatl(Cov_chol));
      Cov_inv = Cov_chol_inv.t() * Cov_chol_inv;
      // Gaussian processes parameters updates
      params = params_aux;
      mgp_sd = mgp_sd_aux;
      mgp_phi = mgp_phi_aux;
      theta_prior_Sigma_chol_inv = theta_prior_Sigma_chol_inv_aux;
      theta_prior_Sigma_inv =
        theta_prior_Sigma_chol_inv.t() * theta_prior_Sigma_chol_inv;
    }

    // Update scaling parameter
    gamma = C / pow(i+1, alpha);
    logscale += gamma * (std::min(exp(accept), 1.0) - target);
    // Rcpp::Rcout << " scale: " << exp(logscale) << std::endl;

    // Update covariance of parameters with stochastic approximation
    arma::vec params_center = params - params_mean;
    params_mean += gamma * params_center;
    params_cov += gamma * (params_center * params_center.t() - params_cov);
    // Rcpp::Rcout << params_cov << std::endl;


    // Save samples
    theta_samples.col(i) = theta;
    c_samples.col(i) = c;
    a_samples.col(i) = arma::vectorise(LA);
    z_samples.col(i) = z;
    mgp_sd_samples.col(i) = mgp_sd;
    mgp_phi_samples.col(i) = mgp_phi;
    corr_samples.col(i) = trimatl2vec(Corr_chol * Corr_chol.t(), false);

  }

    Rcpp::Rcout << "acceptance " << exp(accept) << std::endl;

  return Rcpp::List::create(
      Rcpp::Named("theta") = theta_samples.t(),
      Rcpp::Named("c") = c_samples.t(),
      Rcpp::Named("a") = a_samples.t(),
      Rcpp::Named("z") = z_samples.t(),
      Rcpp::Named("mgp_sd") = mgp_sd_samples.t(),
      Rcpp::Named("mgp_phi") = mgp_phi_samples.t()
      // ,
      // Rcpp::Named("mgp_sd") = mgp_sd,
      // Rcpp::Named("a_Sigma_chol") = a_Sigma_chol,
      // Rcpp::Named("a_mu") = a_mu,
      // Rcpp::Named("a2") = a
      );
  // return Rcpp::List::create(
  //     Rcpp::Named("z_mu") = z_mu
  //     // Rcpp::Named("LA") = LA
  //     );
}

