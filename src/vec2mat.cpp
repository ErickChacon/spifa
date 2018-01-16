#include <RcppArmadillo.h>

//' @export
// [[Rcpp::export]]
arma::mat vec2mat(arma::vec x, int nrow, int ncol) {
  arma::mat X = arma::mat(x);
  X.reshape(nrow, ncol);
  return X;
}

//' @export
// [[Rcpp::export]]
arma::mat vec2matt(arma::vec x, int nrow, int ncol) {
  arma::mat X = arma::mat(x);
  X.reshape(nrow, ncol);
  return X.t();
}

//' @export
// [[Rcpp::export]]
arma::mat vec2ma(arma::vec x, int nrow, int ncol) {
  arma::mat X = arma::mat(x);
  X.reshape(nrow, ncol);
  return X.t();
}


// //' @export
// // [[Rcpp::export]]
// arma::mat subset_cpp(arma::vec x, int nrow, int ncol) {
//   arma::mat X = arma::mat(x);
//   X.reshape(nrow, ncol);
//   return X.t();
// }
//
// //' @export
// // [[Rcpp::export]]
// Rcpp::List subset_cpp(arma::mat X, arma::vec y) {
//   arma::vec plop;
//   plop = X.col(1);
//   return Rcpp::List::create(
//       Rcpp::Named("plop") = plop
//       );
// }
// undefined symbol: _spmirt_subset
