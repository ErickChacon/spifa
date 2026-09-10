
#include <RcppArmadillo.h>
#include "arma-mat.h"
#include "pdf.h"
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppArmadillo)]]

double deviance(arma::vec y, arma::vec c, arma::vec a, arma::vec theta, arma::mat L,
    int n, int q, int m)
{

  // constans
  arma::vec ones_n = arma::ones(n);

  // Compute deviance
  arma::mat Theta = vec2mat(theta, n, m);
  arma::mat LA = L % vec2mat(a, q, m);
  arma::vec linear_pred = arma::kron(c, ones_n) + arma::vectorise(Theta * LA.t());
  arma::vec prob = arma::normcdf(linear_pred);
  arma::vec aux = y % log(prob) + (1 - y) % log(1-prob);
  double deviance = -2 * accu(aux.elem(find_finite(aux)));

  return deviance;
}

// Not exported to R users: called internally by dic.spifa() (R/spifa-methods.R).
// [[Rcpp::export]]
Rcpp::List dic_cpp(arma::vec y, arma::mat c, arma::mat a, arma::mat theta,
    int n, int q, int m, arma::mat L) {

  // number of samples
  int nsamples = c.n_rows;

  // transpose samples for easy access through columns
  c = c.t();
  a = a.t();
  theta = theta.t();

  // average of the parameters by row
  arma::vec c_mean = arma::mean(c, 1);
  arma::vec a_mean = arma::mean(a, 1);
  arma::vec theta_mean = arma::mean(theta, 1);

  // deviance of the posterior average
  double deviance_of_average = deviance(y, c_mean, a_mean, theta_mean, L, n, q, m);

  // posterior average of the deviance
  double mean_deviance = 0;
  for (int i = 0; i < nsamples; ++i) {
    mean_deviance += deviance(y, c.col(i), a.col(i), theta.col(i), L, n, q, m);
  }
  mean_deviance /= nsamples;

  // effective number of parameters and dic
  double p_eff = mean_deviance - deviance_of_average;
  double dic = mean_deviance + p_eff;

  // output list
  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("mean_deviance") = mean_deviance,
      Rcpp::Named("p_eff") = p_eff,
      Rcpp::Named("dic") = dic
      );

  return output;
}


