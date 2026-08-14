//
// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
//
// //' @export
// // [[Rcpp::export]]
// double dinvwish(double v, arma::mat X, arma::mat S) {
//   const int p = X.n_cols;
//   const double pi = M_PI;
//   arma::mat S_chol = arma::chol(S, "lower");
//   double loglike = v * arma::accu(log(S_chol.diag()));
//   arma::mat X_chol = arma::trimatl(arma::chol(X, "lower"));
//   loglike -= (v + p + 1) * arma::accu(log(X_chol.diag()));
//   arma::mat X_chol_inv = X_chol.i();
//   arma::mat X_inv = X_chol_inv.t() * X_chol_inv;
//   loglike -= 0.5 * arma::accu(S % X_inv);
//   loglike -=  0.5 * v * p * log(2);
//   double loggamma = arma::accu(arma::lgamma((v + 1 - arma::linspace(1, p, p)) / 2));
//   loglike -= loggamma + 0.25 * p * (p - 1) * log(pi);
//   return exp(loglike);
// }
//
