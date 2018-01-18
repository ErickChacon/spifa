#include <RcppArmadillo.h>

#ifndef VEC2MAT_H
#define VEC2MAT_H

arma::mat vec2mat(arma::vec x, int nrow, int ncol);
arma::mat vec2matt(arma::vec x, int nrow, int ncol);
arma::mat vec2ma(arma::vec x, int nrow, int ncol);
arma::vec vecsub(arma::vec x, int first_index, int n_length);

#endif

