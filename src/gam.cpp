
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
arma::vec rimvnorm_Q_eig(arma::mat Q) {

  double tol = 1.5e-10;
  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, Q);

  arma::uvec keep = find(abs(eigval) > tol);
  Rcpp::Rcout << "Number of eigenvalues equals to zero: " << eigval.n_elem - keep.n_elem << std::endl;

  eigval = eigval.elem(keep);
  eigvec = eigvec.cols(keep);

  arma::vec x = pow(eigval, - 0.5) % arma::randn(eigval.n_elem);
  x = eigvec * x;

  return x;
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
arma::vec rmvnorm_Q(arma::vec mu, arma::mat Q) {

  const int n = mu.n_elem;

  // unrestricted simulation
  arma::mat Q_chol = arma::chol(Q, "lower");
  arma::vec x = mu + arma::solve(arma::trimatu(Q_chol.t()), arma::randn(n));

  return x;
}



//' @export
// [[Rcpp::export]]
double plop() {
  double x = R::rgamma(0.001, 1 / 0.001);
  return x;
}

//' @export
// [[Rcpp::export]]
Rcpp::List gamcpp(arma::vec y, arma::mat X, arma::mat D, double sigma2, double tau2,
    int niter) {

  const int n = X.n_rows;
  const int p = X.n_cols;
  arma::mat XtX = X.t() * X;
  arma::mat DtD = D.t() * D;
  // const int niter = 1000;

  const double sigma2_a_prior = 0.001;
  const double sigma2_b_prior = 0.001;

  const double tau2_a_prior = 0.001;
  const double tau2_b_prior = 0.001;

  arma::mat beta_samples(p, niter);
  arma::mat sigma2_samples(1, niter);
  arma::mat tau2_samples(1, niter);

  for (int i = 0; i < niter; ++i) {

    // sample beta
    arma::mat beta_Q = XtX / sigma2 + DtD / tau2;
    arma::mat beta_Q_chol = arma::chol(beta_Q, "lower");
    arma::mat beta_mean = solve_sympd_chol(beta_Q_chol, X.t() * y / sigma2);
    arma::vec beta = beta_mean;
    beta += arma::solve(arma::trimatu(beta_Q_chol.t()), arma::randn(p));
    beta_samples.col(i) = beta;

    // sample sigma2
    double sigma2_a = sigma2_a_prior + n / 2;
    double sigma2_b = sigma2_b_prior + arma::accu(square(y - X * beta)) / 2;
    sigma2 = 1 / R::rgamma(sigma2_a, 1 / sigma2_b);
    sigma2_samples.col(i) = sigma2;

    // sample sigma2
    double tau2_a = tau2_a_prior + p / 2 - 1;
    double tau2_b = tau2_b_prior + arma::accu(square(D * beta)) / 2;
    tau2 = 1 / R::rgamma(tau2_a, 1 / tau2_b);
    tau2_samples.col(i) = tau2;

  }

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("beta") = beta_samples.t(),
      Rcpp::Named("sigma2_samples") = sigma2_samples.t(),
      Rcpp::Named("tau2_samples") = tau2_samples.t()
      );

  return output;
}

