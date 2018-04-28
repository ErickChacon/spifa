
#include <RcppArmadillo.h>
#include <iostream>
#include "Power.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

// [[Rcpp::export]]
double classtest(int k, double x)
{
  // Using the class without pointers
  Power powObj(k);
  double Power = powObj.xToTheK(x);
  return Power;
}

