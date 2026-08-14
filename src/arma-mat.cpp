#include <RcppArmadillo.h>

// Internal helpers used by the Ifa sampler (src/ifa.cpp). Not exported to R:
// they are generic matrix-reshaping utilities, not part of the package's
// public API.

arma::mat vec2mat(arma::vec x, int nrow, int ncol) {
  arma::mat X = arma::mat(x);
  X.reshape(nrow, ncol);
  return X;
}

arma::mat vec2matt(arma::vec x, int nrow, int ncol) {
  arma::mat X = arma::mat(x);
  X.reshape(nrow, ncol);
  return X.t();
}

// Applies the m x g restriction matrix T to the block-diagonal covariance of
// g independent Gaussian processes, producing the m x m (times n individuals)
// covariance of the restricted multivariate Gaussian process (see paper
// Section "Allowing further flexibility on the multivariate spatial
// structure").
arma::mat TST(arma::mat mgp_Sigma, arma::mat T) {
  const int m = T.n_rows;              // number of factors
  const int g = T.n_cols;              // number of Gaussian processes
  const int ng = mgp_Sigma.n_rows;     // individuals times Gaussian processes
  const int n = ng/g;                  // number of individuals
  const int nm = n*m;                  // individuals times number of factors
  arma::mat output = arma::zeros(nm, nm);

  // upper triangular
  for (int i = 0; i < (m-1); ++i) {
    for (int j = i + 1; j < m; ++j) {
      for (int k = 0; k < g; ++k) {
        if (T(i,k) != 0 && T(j,k) != 0) {
          output.submat(i*n, j*n, arma::size(n, n)) +=
            T(i,k) * T(j,k) * mgp_Sigma.submat(k*n, k*n, arma::size(n,n));
        }
      }
    }
  }

  // diagonal
  for (int i = 0; i < m; ++i) {
    for (int k = 0; k < g; ++k) {
      if (T(i,k) != 0) {
        output.submat(i*n, i*n, arma::size(n,n)) +=
          pow(T(i,k), 2) * mgp_Sigma.submat(k*n, k*n, arma::size(n,n));
      }
    }
  }

  // lower triangular
  for (int i = 1; i < m; ++i) {
    for (int j = 0; j < i; ++j) {
      output.submat(i*n, j*n, arma::size(n, n)) +=
        output.submat(j*n, i*n, arma::size(n,n));
    }
  }

  return output;
}
