
#include <RcppArmadillo.h>

// Internal link-function helpers used by src/ifa-dic.cpp. Not exported to R:
// not part of the package's public API.

double logit(double p) {
  return log(p/(1-p));
}

double logistic(double x) {
  return 1.0 / (1.0 + exp(-x));
}
