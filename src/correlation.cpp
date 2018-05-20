
#include <RcppArmadillo.h>
#include "arma-mat.h"
// [[Rcpp::depends(RcppArmadillo)]]
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]

//' @export
// [[Rcpp::export]]
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


//' @export
// [[Rcpp::export]]
arma::mat vec2trimatl_old(arma::vec x, int K, bool diag = true) {
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
arma::mat vec2trimatl_test(arma::vec x, int K, bool diag = true) {
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
  // omp_set_num_threads(3);
  x = tanh(x);
  arma::mat L = vec2trimatl(x, K, false);
  arma::mat L_chol(K, K, arma::fill::zeros);

  L_chol.col(0) = L.col(0);
  L_chol(0,0) = 1.0;

  // #pragma omp parallel for
  //
  // #pragma omp parallel for num_threads(4)
  for (int i = (K-1); i > 1; --i) {
  // for (int i = 2; i < K; ++i) {
    // Rcpp::Rcout << i << std::endl;
    for (int j = 1; j < i; ++j) {
      L_chol(i,j) = L(i,j) * sqrt(1 - accu(square(L_chol.submat(i,0, i,j-1))));
    }
  }
  for (int i = 1; i < K; ++i) {
    L_chol(i,i) = sqrt(1 - accu(square(L_chol.submat(i,0, i,i-1))));
  }

  return L_chol;
}

// //' @export
// // [[Rcpp::export]]
// double sum_cpp(arma::vec x)
// {
//   omp_set_num_threads(3);
//   double out = 0;
//   #pragma omp parallel for
//   for (int i = 0; i < x.n_elem; ++i) {
//     // Rcpp::Rcout << i << std::endl;
//     out += x(i);
//   }
//   return out;
// }

//' @export
// [[Rcpp::export]]
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

