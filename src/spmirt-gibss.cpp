
// LOAD LIBRARIES --------------------------------------------------------------

#include <RcppArmadillo.h>
#include <RcppTN.h>
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppArmadillo)]]

// FUNCTION PROTOTYPES ---------------------------------------------------------

arma::vec test(arma::vec y);
Rcpp::List ifa_gibbs(Rcpp::NumericVector y, int n, int q, int N, int m);
arma::mat vec2mat(arma::vec x, int nrow, int ncol);
arma::mat vec2matt(arma::vec x, int nrow, int ncol);

// SIMPLE GIBBS SAMPLING -------------------------------------------------------

//' @export
// [[Rcpp::export]]
Rcpp::List ifa_gibbs(Rcpp::NumericVector y, int n, int q, int N, int m = 1) {

  // arma::sp_mat I_q = arma::speye<arma::sp_mat>(5,5);
  // arma::vec ones_n = arma::ones<arma::vec>(n);
  arma::vec ones_n(n, arma::fill::ones);
  arma::mat eye_q = arma::eye<arma::mat>(q,q);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);

  // Rcpp::NumericVector lows = low_thresh[y];
  Rcpp::NumericVector lows = low_thresh[y];
  Rcpp::NumericVector highs = high_thresh[y];
  // Rcpp::NumericVector a = rtn(0.0, 1.0, lows, highs);

  // Initializing c, a, z

  arma::vec c(q, arma::fill::randn);
  arma::vec mu_c(q);
  arma::mat x_c = arma::kron(eye_q, ones_n);
  arma::mat c_mat;

  arma::vec a(q * m, arma::fill::randn);
  arma::mat Sigma_a(q*m, q*m);
  arma::vec mu_a(q*m);
  arma::mat a_mat;

  arma::vec z = arma::zeros<arma::vec>(y.size());
  arma::vec mu_z(y.size());
  Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itmu_z;
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  arma::mat z_mat;

  arma::vec theta(n*m);
  arma::mat Sigma_theta(n*m, n*m);
  arma::vec mu_theta(n*m);
  arma::mat theta_mat;

  for (int i = 0; i < N; ++i) {

    arma::mat A = vec2matt(a, m, q);
    Sigma_theta = arma::kron(A.t() * A, eye_n) +  arma::eye<arma::mat>(n*m, n*m);
    Sigma_theta = arma::inv(Sigma_theta);
    mu_theta = Sigma_theta * arma::kron(A.t(), eye_n) * (z - x_c * c);
    theta = mu_theta + arma::chol(Sigma_theta, "lower") * arma::randn<arma::vec>(n*m);

    mu_c = x_c.t() * (z - arma::kron(A, eye_n) * theta) / (n + 1);
    c = mu_c + arma::randn<arma::vec>(q) / sqrt(n + 1);

    arma::mat Theta = vec2mat(theta, n, m);
    Sigma_a = arma::kron(eye_q, Theta.t() * Theta) +  arma::eye<arma::mat>(q*m, q*m);
    Sigma_a = arma::inv(Sigma_a);
    mu_a = Sigma_a * arma::kron(eye_q, Theta.t()) * (z - x_c * c);
    a = mu_a + arma::chol(Sigma_a, "lower") * arma::randn<arma::vec>(q*m);


    mu_z = x_c * c + arma::kron(eye_q, Theta) * a;
    ity = y.begin();
    itmu_z = mu_z.begin();
    for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
      *it = RcppTN::rtn1(*itmu_z, 1.0, low_thresh[*ity], high_thresh[*ity]);
      ity++;
      itmu_z++;
    }

    theta_mat.insert_cols(i, theta);
    c_mat.insert_cols(i, c);
    a_mat.insert_cols(i, a);
    z_mat.insert_cols(i, z);

  }

  return Rcpp::List::create(
      // Rcpp::Named("ones") = ones_n,
      // Rcpp::Named("eye") = eye_q,
      // Rcpp::Named("x_c") = x_c,
      // Rcpp::Named("ylog") = y > 0,
      // Rcpp::Named("highs") = highs,
      // Rcpp::Named("lows") = lows,
      // Rcpp::Named("y") = test(y),
      // Rcpp::Named("z") = z,
      // Rcpp::Named("mu_z") = mu_z,
      // Rcpp::Named("a") = a,
      // Rcpp::Named("A") = A,
      // Rcpp::Named("theta") = theta,
      Rcpp::Named("theta") = theta_mat.t(),
      Rcpp::Named("c") = c_mat.t(),
      Rcpp::Named("a") = a_mat.t(),
      Rcpp::Named("z") = z_mat.t()
      );
}

//' @export
// [[Rcpp::export]]
Rcpp::List ifa_gibbs_no(Rcpp::NumericVector y, int n, int q, int N, int m = 1) {

  // arma::sp_mat I_q = arma::speye<arma::sp_mat>(5,5);
  // arma::vec ones_n = arma::ones<arma::vec>(n);
  arma::vec ones_n(n, arma::fill::ones);
  arma::mat eye_q = arma::eye<arma::mat>(q,q);
  arma::mat eye_n = arma::eye<arma::mat>(n,n);
  Rcpp::NumericVector low_thresh = Rcpp::NumericVector::create(R_NegInf, 0);
  Rcpp::NumericVector high_thresh = Rcpp::NumericVector::create(0, R_PosInf);

  // Rcpp::NumericVector lows = low_thresh[y];
  Rcpp::NumericVector lows = low_thresh[y];
  Rcpp::NumericVector highs = high_thresh[y];
  // Rcpp::NumericVector a = rtn(0.0, 1.0, lows, highs);

  // Initializing c, a, z

  arma::vec c(q, arma::fill::zeros);
  // arma::vec c(q, arma::fill::randn);
  // arma::vec mu_c(q);
  arma::mat x_c = arma::kron(ones_n, eye_q);
  arma::mat c_mat;

  arma::vec a(q * m, arma::fill::ones);
  // arma::vec a(q * m, arma::fill::randn);
  // arma::mat Sigma_a(q*m, q*m);
  // arma::vec mu_a(q*m);
  arma::mat a_mat;

  arma::vec z = arma::zeros<arma::vec>(y.size());
  arma::vec mu_z(y.size());
  Rcpp::NumericVector::iterator ity = y.begin();
  Rcpp::NumericVector::iterator itmu_z;
  for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
    *it = RcppTN::rtn1(0.0, 1.0, low_thresh[*ity], high_thresh[*ity]);
    ity++;
  }
  arma::mat z_mat;

  arma::mat Sigma_theta(n*m, n*m);
  arma::vec mu_theta(n*m);
  arma::vec theta(n*m);
  arma::mat theta_mat;

  for (int i = 0; i < N; ++i) {

    arma::mat A = vec2mat(a, q, m);
    Sigma_theta = arma::kron(A.t() * A, eye_n) +  arma::eye<arma::mat>(n*m, n*m);
    Sigma_theta = arma::inv(Sigma_theta);
    mu_theta = Sigma_theta * arma::kron(A.t(), eye_n) * (z - x_c * c);
    theta = mu_theta + arma::chol(Sigma_theta, "lower") * arma::randn<arma::vec>(n*m);

    // mu_c = x_c.t() * (z - arma::kron(A, eye_n) * theta) / (n + 1);
    // c = mu_c + arma::randn<arma::vec>(q) / sqrt(n + 1);

    arma::mat Theta = vec2mat(theta, n, m);
    // Sigma_a = arma::kron(Theta.t() * Theta, eye_q) +  arma::eye<arma::mat>(q*m, q*m);
    // Sigma_a = arma::inv(Sigma_a);
    // mu_a = Sigma_a * arma::kron(Theta.t(), eye_q) * (z - x_c * c);
    // a = mu_a + arma::chol(Sigma_a, "lower") * arma::randn<arma::vec>(q*m);


    mu_z = x_c * c + arma::kron(Theta, eye_q) * a;
    ity = y.begin();
    itmu_z = mu_z.begin();
    for (arma::vec::iterator it = z.begin(); it != z.end() ; ++it) {
      *it = RcppTN::rtn1(*itmu_z, 1.0, low_thresh[*ity], high_thresh[*ity]);
      ity++;
      itmu_z++;
    }

    theta_mat.insert_cols(i, theta);
    c_mat.insert_cols(i, c);
    a_mat.insert_cols(i, a);
    z_mat.insert_cols(i, z);

  }

  return Rcpp::List::create(
      // Rcpp::Named("ones") = ones_n,
      // Rcpp::Named("eye") = eye_q,
      // Rcpp::Named("x_c") = x_c,
      // Rcpp::Named("ylog") = y > 0,
      // Rcpp::Named("highs") = highs,
      // Rcpp::Named("lows") = lows,
      // Rcpp::Named("y") = test(y),
      // Rcpp::Named("z") = z,
      // Rcpp::Named("mu_z") = mu_z,
      // Rcpp::Named("a") = a,
      // Rcpp::Named("A") = A,
      // Rcpp::Named("theta") = theta,
      Rcpp::Named("theta") = theta_mat.t(),
      Rcpp::Named("c") = c_mat.t(),
      Rcpp::Named("a") = a_mat.t(),
      Rcpp::Named("z") = z_mat.t()
      );
}

//' @export
// [[Rcpp::export]]
arma::vec test(arma::vec y) {
  return y;
}

// //' @export
// // [[Rcpp::export]]
// arma::mat a2mat(arma::vec a, int q, int m) {
//   arma::mat A = arma::mat(a);
//   A.reshape(q, m);
//   return A;
// }

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
arma::mat theta2mat(arma::vec a, int q, int m) {
  arma::mat A = arma::mat(a);
  A.reshape(q, m);
  return A;
}









// // [[Rcpp::export]]
// List ifa_gibbs() {
//
//     CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
//     NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
//     List z            = List::create( x, y ) ;
//
//     return z ;
// }

