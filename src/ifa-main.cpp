
#include <RcppArmadillo.h>
#include <RcppTN.h>
#include "arma-mat.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppTN)]]

class Ifa {
private:

  // Metadata
  const int n, q, m, n_corr;
  // Data
  Rcpp::NumericVector y;
  // Parameters to sample
  arma::vec c;                        // difficulty
  arma::vec a;                        // discrimination
  arma::vec theta;                    // latent abilities
  arma::vec z;                        // augmented variable
  // Restrictions
  const arma::mat L;                  // restrictions on discrimination
  // Parameters-dependent objects
  arma::mat LA;                       // restricted discrimation
  // Constant objects usefull for sampling
  const arma::vec ones_n;
  const arma::mat eye_q;
  const arma::mat eye_n;
  const arma::mat eye_m;
  const arma::vec zeros_nm;
  const Rcpp::NumericVector low_thresh ;
  const Rcpp::NumericVector high_thresh;

public:

  Ifa(Rcpp::NumericVector response, int nobs, int nitems, int nfactors,
      arma::mat L_rest);
  void update_theta();
  void update_c();
  void update_a();
  void update_z();
  Rcpp::List sample(int niter);
  arma::vec get();
};

Ifa::Ifa (Rcpp::NumericVector response, int nobs, int nitems, int nfactors,
    arma::mat L_rest):
  y(response), n(nobs), q(nitems), m(nfactors),
  n_corr((m-1) * m / 2),
  ones_n(arma::ones(n)),
  zeros_nm(arma::zeros(n*m)),
  eye_q(arma::eye(q,q)),
  eye_n(arma::eye(n,n)),
  eye_m(arma::eye(m,m)),
  low_thresh(Rcpp::NumericVector::create(R_NegInf, 0)),
  high_thresh(Rcpp::NumericVector::create(0, R_PosInf)),
  L(L_rest)
{

  // Initializing c: difficulty parameters
  c = arma::randn(q) * 0.3;
  // Initializing a: discrimation parameters
  a = arma::ones(q*m) * 0.5;
  LA = L % vec2matt(a, m, q);
  // Initializing z: auxiliary variables
  z = arma::zeros(y.size());
  Rcpp::NumericVector::iterator ity = y.begin();
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  // Initializing theta: latent abilities
  theta = arma::vec(q, m).fill(NA_REAL);

}

void Ifa::update_theta(){
  arma::mat aux_Sigma_inv_chol = arma::chol(LA.t() * LA + eye_m, "lower");
  arma::mat aux_Sigma_chol = arma::inv(trimatl(aux_Sigma_inv_chol)).t();
  arma::mat aux_Sigma = aux_Sigma_chol * aux_Sigma_chol.t();
  arma::mat residual = vec2mat(z - arma::kron(c, ones_n), n, q);
  arma::vec theta_mu = arma::vectorise(residual * LA * aux_Sigma);
  theta = theta_mu + arma::kron(aux_Sigma, eye_n) * arma::randn<arma::vec>(n*m);
}

void Ifa::update_c() {
  arma::mat Theta = vec2mat(theta, n, m);
  arma::mat c_X = arma::kron(eye_q, ones_n);
  arma::vec c_mu = c_X.t() * (z - arma::vectorise(Theta * LA.t())) / (n + 1);
  c = c_mu + arma::randn<arma::vec>(q) / sqrt(n + 1);
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
  arma::vec a_mu = a_Sigma_Lt * arma::kron(eye_q, Theta.t()) *
    (z - arma::kron(c, ones_n));
  a = a_mu + a_Sigma_chol * arma::randn<arma::vec>(q*m);
  arma::mat A = vec2matt(a, m, q);
  LA = L % A;
}

void Ifa::update_z() {
  arma::mat Theta = vec2mat(theta, n, m);
  arma::vec z_mu = arma::kron(c, ones_n) + arma::vectorise(Theta * LA.t());
  Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itz_mu = y.begin();
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

  return Rcpp::List::create(
      Rcpp::Named("theta") = theta_samples.t(),
      Rcpp::Named("c") = c_samples.t(),
      Rcpp::Named("a") = a_samples.t(),
      Rcpp::Named("z") = z_samples.t()
      );
}

arma::vec Ifa::get() {
  return zeros_nm;
}

// [[Rcpp::export]]
Rcpp::List spmirt(Rcpp::NumericVector response, int nobs, int nitems, int nfactors,
    arma::mat L_rest, int niter) {

  Ifa model(response, nobs, nitems, nfactors, L_rest);
  Rcpp::List output = model.sample(niter);
  return output;
}


  // if (n*q != response.size()) {
  //   Rcpp::stop("length of 'response' should be equal to 'nobs * nitems'");
  // }
