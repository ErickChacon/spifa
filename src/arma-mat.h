#include <RcppArmadillo.h>

#ifndef ARMA_MAT_H
#define ARMA_MAT_H

arma::mat vec2mat(arma::vec x, int nrow, int ncol);
arma::mat vec2matt(arma::vec x, int nrow, int ncol);
arma::mat TST(arma::mat mgp_Sigma, arma::mat T);

#endif

