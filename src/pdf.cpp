
#include <RcppArmadillo.h>
#include "correlation.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
double dmvnorm(arma::mat X, arma::mat Mean, arma::mat Sigma, bool logpdf = true) {
  const int q = X.n_rows;
  const int n = X.n_cols;
  const double pi = M_PI;
  arma::mat L = arma::chol(Sigma, "lower");
  arma::mat kern_sq_root = arma::solve(arma::trimatl(L), X - Mean);
  double log_pdf = - 0.5 * n * q * log(2.0 * pi);
  log_pdf -= n * arma::accu(log(L.diag()));
  log_pdf -= 0.5 * arma::accu(square(kern_sq_root));
  return (logpdf) ? log_pdf: exp(log_pdf);
}

//' @export
// [[Rcpp::export]]
double dmvnorm_chol(arma::mat X, arma::mat Mean, arma::mat L, bool logpdf = true) {
  const int q = X.n_rows;
  const int n = X.n_cols;
  const double pi = M_PI;
  arma::mat kern_sq_root = arma::solve(arma::trimatl(L), X - Mean);
  double log_pdf = - 0.5 * n * q * log(2.0 * pi);
  log_pdf -= n * arma::accu(log(L.diag()));
  log_pdf -= 0.5 * arma::accu(square(kern_sq_root));
  return (logpdf) ? log_pdf: exp(log_pdf);
}

//' @export
// [[Rcpp::export]]
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

//' @export
// [[Rcpp::export]]
double dmvnorm_prec(arma::vec x, arma::vec mean, arma::mat sigma_inv) {
  const int n = x.n_elem;
  const double pi = M_PI;
  double val;
  double sign;
  arma::log_det(val, sign, sigma_inv);
  double kern = arma::as_scalar((x - mean).t() * sigma_inv * (x - mean));
  double log_pdf = - 0.5 * n * log(2.0 * pi);
  log_pdf += + 0.5 * val;
  log_pdf += - 0.5 * kern;
  return log_pdf;
}

//' @export
// [[Rcpp::export]]
double dinvwish(double v, arma::mat X, arma::mat S, bool logpdf = true) {
  const int p = X.n_cols;
  const double pi = M_PI;
  arma::mat S_chol = arma::chol(S, "lower");
  double log_pdf = v * arma::accu(log(S_chol.diag()));
  arma::mat X_chol = arma::trimatl(arma::chol(X, "lower"));
  log_pdf -= (v + p + 1) * arma::accu(log(X_chol.diag()));
  arma::mat X_chol_inv = X_chol.i();
  arma::mat X_inv = X_chol_inv.t() * X_chol_inv;
  log_pdf -= 0.5 * arma::accu(S % X_inv);
  log_pdf -=  0.5 * v * p * log(2);
  double loggamma = arma::accu(arma::lgamma((v + 1 - arma::linspace(1, p, p)) / 2));
  log_pdf -= loggamma + 0.25 * p * (p - 1) * log(pi);
  return (logpdf) ? log_pdf: exp(log_pdf);
}

//' @export
// [[Rcpp::export]]
double dlkj_corr(arma::mat R, double eta, bool logpdf = true) {
  const int K = R.n_rows;
  arma::mat L = arma::chol(R);
  double log_pdf = 2 * (eta - 1) * arma::accu(log(L.diag()));
  return (logpdf) ? log_pdf: exp(log_pdf);
}

//' @export
// [[Rcpp::export]]
double dlkj_corr_chol(arma::mat L, double eta, bool logpdf = true) {
  const int K = L.n_rows;
  double log_pdf = 0;
  for (int i = 1; i < K; ++i) {
    log_pdf += (K - i-1 + 2*eta-2) * log(L(i, i));
  }
  return (logpdf) ? log_pdf: exp(log_pdf);
}

//' @export
// [[Rcpp::export]]
double dlkj_corr_free(arma::vec x, int K, double eta, bool logpdf = true) {
  // There is problems when some element on x is equal to zero
  arma::mat L = vec2trimatl(tanh(x), K, false);
  arma::mat L_chol = vec2chol_corr(x, K);
  double log_pdf = dlkj_corr_chol(L_chol, eta, true);
  log_pdf -= 2 * arma::accu(log(cosh(x)));
  log_pdf += arma::accu(log(trimatl2vec(L_chol / L, false)));
  return (logpdf) ? log_pdf: exp(log_pdf);
}

//' @export
// [[Rcpp::export]]
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


//' @export
// [[Rcpp::export]]
arma::mat test1(arma::vec A, arma::mat B) {
  B.each_col() %= A;
  // B.each_col() %= A;
  return B;
}

//' @export
// [[Rcpp::export]]
Rcpp::List testing(arma::mat X, arma::vec y) {
  X(1,1) = 1000;
  y(1) = 1000;
  arma::vec di = log(X.diag());
  double plop = arma::accu(y);
  arma::vec plop2 = square(y);
  arma::mat Sigma_proposal(2,2, arma::fill::zeros);
  Sigma_proposal(0,0) = 0.1;
  Sigma_proposal(1,1) = 0.1;
  // bool a = arma::as_scalar(arma::randu<arma::vec>(1)) > 0.5;
  bool a = R::runif(0,1) > 0.5;
  arma::vec scal = y(1) * y;

  // arma::mat L = arma::trimatl(A);
  //
  return Rcpp::List::create(
      Rcpp::Named("y") = y,
      Rcpp::Named("X") = X,
      Rcpp::Named("di") = di,
      Rcpp::Named("plop") = plop,
      Rcpp::Named("plop2") = plop2,
      Rcpp::Named("Sigma_proposal") = Sigma_proposal,
      Rcpp::Named("a") = a,
      Rcpp::Named("scal") = scal
      );
}


