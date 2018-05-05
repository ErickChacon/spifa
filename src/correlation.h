#include <RcppArmadillo.h>

#ifndef CORRELATION_H
#define CORRELATION_H

arma::mat vec2trimatl(arma::vec x, int K, bool diag = true);
arma::vec trimatl2vec(arma::mat L, bool diag = true);
arma::mat vec2chol_corr(arma::vec x, int K);
Rcpp::List vec2chol_corr2(arma::vec x, int K);
arma::vec chol_corr2vec(arma::mat L_chol);

#endif


