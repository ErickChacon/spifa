#include <RcppArmadillo.h>

#ifndef PDF_H
#define PDF_H

double dmvnorm_cholinv(arma::mat X, arma::mat Mean, arma::mat L_inv,
    bool logpdf = true);
double dlkj_corr_chol(arma::mat L, double eta, bool logpdf = true);
double dlkj_corr_free2(arma::vec x, int K, double eta, bool logpdf = true);

#endif

