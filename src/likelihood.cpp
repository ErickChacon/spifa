
#include <RcppArmadillo.h>

//' @export
// [[Rcpp::export]]
double dmvnorm(arma::vec x, arma::vec mean, arma::mat sigma) {
  const int n = x.n_elem;
  const double pi = M_PI;
  arma::mat L = arma::chol(sigma, "lower");
  arma::mat kern_sq_root = arma::solve(arma::trimatl(L), x - mean);
  double loglike = - 0.5 * n * log(2.0 * pi);
  loglike += - arma::accu(log(L.diag()));
  loglike += - 0.5 * arma::accu(square(kern_sq_root));
  return loglike;
}

//' @export
// [[Rcpp::export]]
double dmvnorm_chol(arma::vec x, arma::vec mean, arma::mat L) {
  const int n = x.n_elem;
  const double pi = M_PI;
  arma::mat kern_sq_root = arma::solve(arma::trimatl(L), x - mean);
  double loglike = - 0.5 * n * log(2.0 * pi);
  loglike += - arma::accu(log(L.diag()));
  loglike += - 0.5 * arma::accu(square(kern_sq_root));
  return loglike;
}

//' @export
// [[Rcpp::export]]
double dmvnorm_cholinv(arma::vec x, arma::vec mean, arma::mat L_inv) {
  const int n = x.n_elem;
  const double pi = M_PI;
  arma::mat kern_sq_root = L_inv * (x - mean);
  double loglike = - 0.5 * n * log(2.0 * pi);
  loglike += arma::accu(log(L_inv.diag()));
  loglike += - 0.5 * arma::accu(square(kern_sq_root));
  return loglike;
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
  double loglike = - 0.5 * n * log(2.0 * pi);
  loglike += + 0.5 * val;
  loglike += - 0.5 * kern;
  return loglike;
}

//' @export
// [[Rcpp::export]]
double dinvwish(double v, arma::mat X, arma::mat S) {
  const int p = X.n_cols;
  const double pi = M_PI;
  arma::mat S_chol = arma::chol(S, "lower");
  double loglike = v * arma::accu(log(S_chol.diag()));
  arma::mat X_chol = arma::trimatl(arma::chol(X, "lower"));
  loglike -= (v + p + 1) * arma::accu(log(X_chol.diag()));
  arma::mat X_chol_inv = X_chol.i();
  arma::mat X_inv = X_chol_inv.t() * X_chol_inv;
  loglike -= 0.5 * arma::accu(S % X_inv);
  loglike -=  0.5 * v * p * log(2);
  double loggamma = arma::accu(arma::lgamma((v + 1 - arma::linspace(1, p, p)) / 2));
  loglike -= loggamma + 0.25 * p * (p - 1) * log(pi);
  return exp(loglike);
}

//' @export
// [[Rcpp::export]]
double dlkj_corr(arma::mat R, double eta) {
  const int K = R.n_rows;
  arma::mat L = arma::chol(R);
  double loglike = 2 * (eta - 1) * arma::accu(log(L.diag()));
  return exp(loglike);
}

//' @export
// [[Rcpp::export]]
double dlkj_corr_chol(arma::mat L, double eta) {
  const int K = L.n_rows;
  double loglike = 0;
  for (int i = 1; i < K; ++i) {
    loglike += (K - i-1 + 2*eta-2) * log(L(i, i));
  }
  return exp(loglike);
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


