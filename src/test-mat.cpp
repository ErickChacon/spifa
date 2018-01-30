#include <RcppArmadillo.h>

//' @export
// [[Rcpp::export]]
arma::mat matmatmat(arma::mat x) {
  return x * x * x;
}

//' @export
// [[Rcpp::export]]
arma::mat matmat(arma::mat x) {
  return x * x;
}

//' @export
// [[Rcpp::export]]
arma::mat ar_chol(arma::mat x) {
  return arma::chol(x);
}


//' @export
// [[Rcpp::export]]
arma::mat mat_inv(arma::mat x) {
  return x.i();
}
//' @export
// [[Rcpp::export]]
arma::mat mat_inv2(arma::mat x) {
  return arma::inv_sympd(x);
}

//' @export
// [[Rcpp::export]]
arma::mat inv_chol(arma::mat x) {
  arma::mat a = arma::solve(trimatl(arma::chol(x, "lower")),
      arma::eye<arma::mat>(size(x)));
  return a.t() * a;
}

// //' @export
// // [[Rcpp::export]]
// arma::mat inv_chol2(arma::mat x) {
//   arma::mat a = arma::solve((arma::chol(x, "lower")),
//       arma::eye<arma::mat>(size(x)));
//   return a.t() * a;
// }


//' @export
// [[Rcpp::export]]
arma::mat inv_solve(arma::mat x, arma::mat eye) {
  return (x.i() + eye).i();
}

//' @export
// [[Rcpp::export]]
arma::mat inv_solve2(arma::mat x, arma::mat eye) {
  return arma::inv_sympd(arma::inv_sympd(x) + eye);
}

//' @export
// [[Rcpp::export]]
arma::mat inv_form(arma::mat x, arma::mat eye) {
  arma::mat inv_aux = arma::inv_sympd(x + eye);
  return x * (eye - inv_aux * x);
}
