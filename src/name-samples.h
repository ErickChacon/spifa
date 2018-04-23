// #include <RcppArmadillo.h>

#ifndef NAME_SAMPLES_H
#define NAME_SAMPLES_H

Rcpp::StringVector name_samples_vec(int n_elem, std::string name);
Rcpp::StringVector name_samples_mat(int nrow, int ncol, std::string name);
Rcpp::StringVector name_samples_lower(int nrow, int ncol, std::string name,
    bool diag = true);

#endif

