
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
Rcpp::List testing(arma::mat X, arma::vec y) {
  X(1,1) = 1000;
  y(1) = 1000;
  arma::vec di = log(X.diag());
  double plop = arma::accu(y);
  arma::vec plop2 = square(y);

  // arma::mat L = arma::trimatl(A);
  //
  return Rcpp::List::create(
      Rcpp::Named("y") = y,
      Rcpp::Named("X") = X,
      Rcpp::Named("di") = di,
      Rcpp::Named("plop") = plop,
      Rcpp::Named("plop2") = plop2
      );
}


