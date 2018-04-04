
#include <RcppArmadillo.h>
#include "arma-mat.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
arma::mat vec2trimatl (arma::vec x, int K, bool diag = true) {
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
      L.submat(i+1,i, K-1,i) = vecsub(x, round(i*K - (i+1)*(i)/2), K-i-1);
    }
  }
  return L;
}

//' @export
// [[Rcpp::export]]
arma::mat vec2corr(arma::vec x, int K) {
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

