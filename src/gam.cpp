
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//

//' @export
// [[Rcpp::export]]
arma::mat solve_sympd(arma::mat A, arma::mat B) {
  // AX = B, where A is symmetric positive definite

  arma::mat A_chol = arma::chol(A, "lower");
  arma::mat V = arma::solve(arma::trimatl(A_chol), B);
  arma::mat X = arma::solve(arma::trimatu(A_chol.t()), V);

  return X;
}

//' @export
// [[Rcpp::export]]
arma::mat solve_sympd_chol(arma::mat A_chol, arma::mat B) {
  // AX = B, where A is symmetric positive definite

  arma::mat V = arma::solve(arma::trimatl(A_chol), B);
  arma::mat X = arma::solve(arma::trimatu(A_chol.t()), V);

  return X;
}

//' @export
// [[Rcpp::export]]
arma::vec rmvnorm_rest_Q(arma::vec mu, arma::mat Q, arma::mat A, arma::vec e) {

  const int n = mu.n_elem;

  // unrestricted simulation
  arma::mat Q_chol = arma::chol(Q, "lower");
  arma::vec x = mu + arma::solve(arma::trimatu(Q_chol.t()), arma::randn(n));

  // transform realization to fullfill contrains
  arma::mat V = solve_sympd_chol(Q_chol, A.t());
  arma::mat W = A * V;
  arma::mat U = solve_sympd(W, V.t());
  arma::vec c = A * x - e;
  x = x - U.t() * c;

  return x;

}

//' @export
// [[Rcpp::export]]
arma::mat gamcpp(arma::vec y, arma::mat X, arma::mat D, double sigma2, double tau2) {

  // const int n = y.n_elem;
  const int p = X.n_cols;
  const int niter = 1000;
  arma::mat XtX = X.t() * X;
  arma::mat DtD = D.t() * D;

  arma::mat beta_samples(p, niter);

  for (int i = 0; i < niter; ++i) {
    arma::mat beta_Q = XtX / sigma2 + DtD / tau2;
    arma::mat beta_Q_chol = arma::chol(beta_Q, "lower");
    arma::mat beta_mean = solve_sympd_chol(beta_Q_chol, X.t() * y) / sigma2;
    arma::vec beta = beta_mean;
    beta += arma::solve(arma::trimatu(beta_Q_chol.t()), arma::randn(p));
    beta_samples.col(i) = beta;
  }

  return beta_samples;
}

