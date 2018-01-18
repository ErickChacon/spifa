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
Rcpp::List subset_cpp(arma::mat X, arma::vec y) {
  arma::vec sub_col = X.col(0);
  arma::vec sub_vec = y.subvec(0, 3);
  double element = arma::as_scalar(y.subvec(0, 0));
  arma::vec z(10, arma::fill::zeros);
  arma::vec::iterator itz = z.begin();
  for (int i = 0; i < 10; ++i) {
    *itz = i;
    itz++;
  }
  return Rcpp::List::create(
      Rcpp::Named("col") = sub_col,
      Rcpp::Named("sub_vec") = sub_vec,
      Rcpp::Named("elem") = element,
      Rcpp::Named("loop_fill") = z
      );
}

