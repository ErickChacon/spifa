
#include <RcppArmadillo.h>

//' @export
// [[Rcpp::export]]
double logit(double p) {
  return log(p/(1-p));
}
// double logit(double p) {
//   return log(p/(1-p));
// }

//' @export
// [[Rcpp::export]]
double logistic(double x) {
  return 1.0 / (1.0 + exp(-x));
}
