
#include <RcppArmadillo.h>
#include "arma-mat.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
arma::mat vec2trimatl(arma::vec x, int K, bool diag = true) {
  // K = (-1 + sqrt(1 + 8 * x.length))/2
  // K = (1 + sqrt(1 + 8 * x.length))/2
  arma::mat L(K, K, arma::fill::zeros);
  // fill lower triangle by column
  if (diag) {
    for (int i = 0; i < K; ++i) {
      L.submat(i,i, K-1,i) = vecsub(x, round(i*K - i*(i-1)/2), K-i);
    }
  } else {
    for (int i = 0; i < K-1; ++i) {
      L.submat(i+1,i, K-1,i) = vecsub(x, round(i*K - (i+1)*i/2), K-i-1);
    }
  }
  return L;
}

//' @export
// [[Rcpp::export]]
arma::vec trimatl2vec(arma::mat L, bool diag = true) {
  int K = L.n_rows;
  int N;

  if (diag) {
    N = round(K * (K + 1) / 2);
  } else {
    N = round((K - 1) * K / 2);
  }

  arma::vec x(N, arma::fill::zeros);
  // fill vector by lower triangle column
  if (diag) {
    for (int i = 0; i < K; ++i) {
      x.subvec(round(i*K - i*(i-1)/2), arma::size(K-i,1)) = L.submat(i,i, K-1,i);
    }
  } else {
    for (int i = 0; i < K-1; ++i) {
      x.subvec(round(i*K-(i+1)*i/2), arma::size(K-i-1, 1)) = L.submat(i+1,i, K-1,i);
    }
  }

  return x;
}


//' @export
// [[Rcpp::export]]
arma::mat vec2chol_corr(arma::vec x, int K) {
  x = tanh(x);
  arma::mat L = vec2trimatl(x, K, false);
  arma::mat L_chol(K, K, arma::fill::zeros);

  L_chol.col(0) = L.col(0);
  L_chol(0,0) = 1.0;

  for (int i = 2; i < K; ++i) {
    for (int j = 1; j < i; ++j) {
      L_chol(i,j) = L(i,j) * sqrt(1 - accu(square(L_chol.submat(i,0, i,j-1))));
    }
  }
  for (int i = 1; i < K; ++i) {
    L_chol(i,i) = sqrt(1 - accu(square(L_chol.submat(i,0, i,i-1))));
  }

  return L_chol;
}


//' @export
// [[Rcpp::export]]
arma::vec chol_corr2vec(arma::mat L_chol) {
  int K = L_chol.n_rows;
  arma::mat L(K, K, arma::fill::zeros);

  L.col(0) = L_chol.col(0);
  L(0,0) = 0;

  for (int i = 2; i < K; ++i) {
    for (int j = 1; j < i; ++j) {
      L(i,j) = L_chol(i,j) / sqrt(1 - accu(square(L_chol.submat(i,0, i,j-1))));
    }
  }

  arma::vec x = trimatl2vec(L, false);
  x = atanh(x);

  return x;
}
