
#include <RcppArmadillo.h>
#include "ifa.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double spmirt(double y) {
  Ifa obj(2);
  return obj.prior(y);
}


