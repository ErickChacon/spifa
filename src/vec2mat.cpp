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
arma::vec vecsub(arma::vec x, int first_index, int n_length) {
  arma::vec sub_x = x.subvec(first_index, first_index + n_length - 1);
  return sub_x;
}

//' @export
// [[Rcpp::export]]
double vecsub1(arma::vec x, int index) {
  double out = arma::as_scalar(x.subvec(index, index));
  return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::List subset_cpp(arma::mat X, arma::vec y) {
  arma::vec sub_col = X.col(0);
  arma::vec sub_vec = y.subvec(0, 3);
  // arma::vec vec_scal = y - y.subvec(0, 0);
  double vec_scal = vecsub1(y, 1);
  // double size_vec = size(y);
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
      Rcpp::Named("vec_scal") = vec_scal,
      Rcpp::Named("elem") = element,
      Rcpp::Named("loop_fill") = z
      );
}

