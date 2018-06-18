
#include <RcppArmadillo.h>
#include "arma-mat.h"
#include "pdf.h"
#include "links.h"
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
double dic_cpp(arma::vec y, arma::mat c, arma::mat a, arma::mat theta,
    int n, int q, int m, arma::mat L) {

  // number of samples
  int nsamples = c.n_rows;

  // transpose samples for easy access through columns
  c = c.t();
  a = a.t();
  theta = theta.t();
  arma::vec ones_n = arma::ones(n);

  // average of the parameters by row
  arma::vec c_mean = arma::mean(c, 1);
  arma::vec a_mean = arma::mean(a, 1);
  arma::vec theta_mean = arma::mean(theta, 1);

  // deviance of the posterior average
  arma::mat Theta_mean = vec2mat(theta_mean, n, m);
  arma::mat LA_mean = L % vec2matt(a_mean, m, q);
  arma::vec eta_mean = arma::kron(c_mean, ones_n) +
    arma::vectorise(Theta_mean * LA_mean.t());
  arma::vec log_prob = log(arma::normcdf(eta_mean));
  log_prob = y % log_prob + (1-y) % log_prob;
  double deviance_of_average = accu(log_prob.elem(find_finite(log_prob)));
  double dic = - deviance_of_average;

  // posterior average of the deviance
  double average_of_deviance = 0;
  for (int i = 0; i < nsamples; ++i) {
    arma::mat Theta = vec2mat(theta.col(i), n, m);
    arma::mat LA = L % vec2matt(a.col(i), m, q);
    arma::vec eta = arma::kron(c.col(i), ones_n) + arma::vectorise(Theta * LA.t());
    arma::vec aux = log(arma::normcdf(eta));
    aux = y % aux + (1-y) % aux;
    average_of_deviance += accu(aux.elem(find_finite(aux)));
  }
  average_of_deviance /= nsamples;
  dic += 2 * average_of_deviance;

  return dic;
}


