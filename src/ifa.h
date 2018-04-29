
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef IFA_H
#define IFA_H

class Ifa {
private:
  int n, q, m;
  arma::vec A;

public:
  Ifa(int n);
  double prior(double y);
  // virtual ~Ifa();
};

#endif /* IFA_H */
