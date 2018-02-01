
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


