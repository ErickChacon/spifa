
#include <RcppArmadillo.h>
#include "arma-mat.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal helpers used by the Ifa sampler (src/ifa.cpp and src/pdf.cpp) to
// move between the unconstrained parameterisation of the residual
// correlation matrix (a real vector, used in the adaptive Metropolis-Hastings
// proposal) and its Cholesky/matrix forms. Not exported to R: they are not
// part of the package's public API.

arma::mat vec2trimatl(arma::vec x, int K, bool diag = true) {
  arma::mat L(K, K, arma::fill::ones);
  if (diag) {
    // fill lower triangle, with diagonal, by column
    L = arma::trimatl(L);
    L.elem(find(trimatl(L))) = x;
  } else {
    // fill lower triangle, without diagonal, by column
    L = arma::trimatl(L, -1);
    L.elem(find(trimatl(L))) = x;
  }
  return L;
}

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

arma::mat vec2chol_corr(arma::vec x, int K) {
  x = tanh(x);
  arma::mat L = vec2trimatl(x, K, false);
  arma::mat L_chol(K, K, arma::fill::zeros);

  L_chol.col(0) = L.col(0);
  L_chol(0,0) = 1.0;

  for (int i = (K-1); i > 1; --i) {
    for (int j = 1; j < i; ++j) {
      L_chol(i,j) = L(i,j) * sqrt(1 - accu(square(L_chol.submat(i,0, i,j-1))));
    }
  }
  for (int i = 1; i < K; ++i) {
    L_chol(i,i) = sqrt(1 - accu(square(L_chol.submat(i,0, i,i-1))));
  }

  return L_chol;
}

Rcpp::List vec2chol_corr2(arma::vec x, int K) {
  x = tanh(x);
  arma::mat L = vec2trimatl(x, K, false);
  arma::mat L_chol(K, K, arma::fill::zeros);
  arma::mat L_grad(K, K, arma::fill::zeros);
  L_grad.col(0) += 1;
  L_grad(0,0) = 0;

  L_chol.col(0) = L.col(0);
  L_chol(0,0) = 1.0;

  for (int i = 2; i < K; ++i) {
    for (int j = 1; j < i; ++j) {
      L_grad(i,j) = sqrt(1 - accu(square(L_chol.submat(i,0, i,j-1))));
      L_chol(i,j) = L(i,j) * L_grad(i,j);
    }
  }
  for (int i = 1; i < K; ++i) {
    L_chol(i,i) = sqrt(1 - accu(square(L_chol.submat(i,0, i,i-1))));
  }

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("L_chol") = L_chol,
      Rcpp::Named("L_grad") = L_grad
      );

  return output;
}

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
