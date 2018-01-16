#include <RcppArmadillo.h>

#ifndef VEC2MAT_H
#define VEC2MAT_H

arma::mat vec2mat(arma::vec x, int nrow, int ncol);
arma::mat vec2matt(arma::vec x, int nrow, int ncol);
arma::mat vec2ma(arma::vec x, int nrow, int ncol);
// arma::mat subset_cpp(arma::vec x, int nrow, int ncol);

#endif

