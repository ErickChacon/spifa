#include <RcppArmadillo.h>

#ifndef PDF_H
#define PDF_H

double dmvnorm(arma::vec x, arma::vec mean, arma::mat sigma);
double dmvnorm_chol(arma::vec x, arma::vec mean, arma::mat sigma_chol);
double dmvnorm_cholinv(arma::vec x, arma::vec mean, arma::mat L_inv);
double dinvwish(double v, arma::mat X, arma::mat S);
double dlkj_corr(arma::mat R, double eta, bool loglik = false);
double dlkj_corr_chol(arma::mat L, double eta, bool loglik = false);
double dlkj_corr_free(arma::vec x, int K, double eta, bool loglik = false);

#endif

