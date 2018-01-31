#include <RcppArmadillo.h>

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

double dmvnorm(arma::vec x, arma::vec mean, arma::mat sigma);
double dmvnorm_chol(arma::vec x, arma::vec mean, arma::mat sigma_chol);

#endif

