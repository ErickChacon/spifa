#include <RcppArmadillo.h>

#ifndef ARMA_MAT_H
#define ARMA_MAT_H

arma::mat vec2trimatl(arma::vec x, int K, bool diag = true);
arma::vec trimatl2vec(arma::mat L, bool diag = true);
arma::mat vec2chol_corr(arma::vec x, int K);
arma::vec chol_corr2vec(arma::mat L_chol);

#endif


