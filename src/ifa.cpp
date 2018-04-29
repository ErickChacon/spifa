
#include <RcppArmadillo.h>
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

Ifa::Ifa (int N) {
  n = N;
}

double Ifa::prior(double y) {
  return pow(y, n);
};

