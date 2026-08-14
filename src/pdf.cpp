
#include <RcppArmadillo.h>
#include "correlation.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal density functions used by the Ifa sampler (src/ifa.cpp) for its
// Gibbs/Metropolis-Hastings updates (multivariate normal, inverse-Wishart,
// and LKJ correlation priors). Not exported to R: they duplicate
// functionality available in packages such as mvtnorm, and are not part of
// the package's public API.

double dmvnorm_cholinv(arma::mat X, arma::mat Mean, arma::mat L_inv,
    bool logpdf = true) {
  const int q = X.n_rows;
  const int n = X.n_cols;
  const double pi = M_PI;
  arma::mat kern_sq_root = L_inv * (X - Mean);
  double log_pdf = - 0.5 * n * q * log(2.0 * pi);
  log_pdf += n * arma::accu(log(L_inv.diag()));
  log_pdf -= 0.5 * arma::accu(square(kern_sq_root));
  return (logpdf) ? log_pdf: exp(log_pdf);
}

double dlkj_corr_chol(arma::mat L, double eta, bool logpdf = true) {
  const int K = L.n_rows;
  double log_pdf = 0;
  for (int i = 1; i < K; ++i) {
    log_pdf += (K - i-1 + 2*eta-2) * log(L(i, i));
  }
  return (logpdf) ? log_pdf: exp(log_pdf);
}

double dlkj_corr_free2(arma::vec x, int K, double eta, bool logpdf = true) {
  arma::mat L = vec2trimatl(tanh(x), K, false);
  Rcpp::List aux = vec2chol_corr2(x, K);
  arma::mat L_chol = aux["L_chol"];
  arma::mat L_grad = aux["L_grad"];
  double log_pdf = dlkj_corr_chol(L_chol, eta, true);
  log_pdf -= 2 * arma::accu(log(cosh(x)));
  log_pdf += arma::accu(log(trimatl2vec(L_grad, false)));
  return (logpdf) ? log_pdf: exp(log_pdf);
}
