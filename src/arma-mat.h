#include <RcppArmadillo.h>

#ifndef ARMA_MAT_H
#define ARMA_MAT_H

arma::mat vec2mat(arma::vec x, int nrow, int ncol);
arma::mat vec2matt(arma::vec x, int nrow, int ncol);
arma::mat vec2ma(arma::vec x, int nrow, int ncol);
arma::vec vecsub(arma::vec x, int first_index, int n_length);
double vecsub1(arma::vec x, int index);
double matsub1(arma::mat x, int index_row, int index_col);
arma::mat TST(arma::mat mgp_Sigma, arma::mat T);

#endif

