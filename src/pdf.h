#include <RcppArmadillo.h>

#ifndef PDF_H
#define PDF_H

double dmvnorm(arma::mat X, arma::mat Mean, arma::mat Sigma, bool logpdf = true);
double dmvnorm_chol(arma::mat X, arma::mat Mean, arma::mat L, bool logpdf = true);
double dmvnorm_cholinv(arma::mat X, arma::mat Mean, arma::mat L_inv,
    bool logpdf = true);
double dinvwish(double v, arma::mat X, arma::mat S, bool logpdf = true);
double dlkj_corr(arma::mat R, double eta, bool logpdf = true);
double dlkj_corr_chol(arma::mat L, double eta, bool logpdf = true);
double dlkj_corr_free(arma::vec x, int K, double eta, bool logpdf = true);
double dlkj_corr_free2(arma::vec x, int K, double eta, bool logpdf = true);

#endif

