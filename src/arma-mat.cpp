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
double matsub1(arma::mat x, int index_row, int index_col) {
  double out = arma::as_scalar(x.submat(index_row,index_col, index_row,index_col));
  return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::List subset_cpp(arma::mat X, arma::vec y) {
  // X(arma::span(0, 0), arma::span(0, 0)) = 5.0;
  X.submat(0, 0, 0, 0) = 5.1;
  X.submat(0, 1, 0, 3) = y.subvec(0,2).t();
  X.diag().ones();
  // X.submat(1, 0, 3, 0) = y.subvec(0,2);
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
  arma::mat A(9,5,arma::fill::randn);
  A.diag() = abs(A.diag());
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 5; ++j) {
      if (j > i) {
        A.submat(i,j,i,j) = 0;
      }
    }
  }
  // arma::mat L = arma::trimatl(A);
  return Rcpp::List::create(
      Rcpp::Named("col") = sub_col,
      Rcpp::Named("sub_vec") = sub_vec,
      Rcpp::Named("vec_scal") = vec_scal,
      Rcpp::Named("elem") = element,
      Rcpp::Named("loop_fill") = z,
      Rcpp::Named("X") = X,
      Rcpp::Named("Lower") = A
      );
}

