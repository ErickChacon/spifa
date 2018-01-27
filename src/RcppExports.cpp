// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// vec2mat
arma::mat vec2mat(arma::vec x, int nrow, int ncol);
RcppExport SEXP _spmirt_vec2mat(SEXP xSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(vec2mat(x, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// vec2matt
arma::mat vec2matt(arma::vec x, int nrow, int ncol);
RcppExport SEXP _spmirt_vec2matt(SEXP xSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(vec2matt(x, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// vecsub
arma::vec vecsub(arma::vec x, int first_index, int n_length);
RcppExport SEXP _spmirt_vecsub(SEXP xSEXP, SEXP first_indexSEXP, SEXP n_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type first_index(first_indexSEXP);
    Rcpp::traits::input_parameter< int >::type n_length(n_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(vecsub(x, first_index, n_length));
    return rcpp_result_gen;
END_RCPP
}
// vecsub1
double vecsub1(arma::vec x, int index);
RcppExport SEXP _spmirt_vecsub1(SEXP xSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(vecsub1(x, index));
    return rcpp_result_gen;
END_RCPP
}
// matsub1
double matsub1(arma::mat x, int index_row, int index_col);
RcppExport SEXP _spmirt_matsub1(SEXP xSEXP, SEXP index_rowSEXP, SEXP index_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type index_row(index_rowSEXP);
    Rcpp::traits::input_parameter< int >::type index_col(index_colSEXP);
    rcpp_result_gen = Rcpp::wrap(matsub1(x, index_row, index_col));
    return rcpp_result_gen;
END_RCPP
}
// subset_cpp
Rcpp::List subset_cpp(arma::mat X, arma::vec y);
RcppExport SEXP _spmirt_subset_cpp(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(subset_cpp(X, y));
    return rcpp_result_gen;
END_RCPP
}
// probit_gp
Rcpp::List probit_gp(Rcpp::NumericVector y, arma::mat dist, double tau2, double phi, int iter);
RcppExport SEXP _spmirt_probit_gp(SEXP ySEXP, SEXP distSEXP, SEXP tau2SEXP, SEXP phiSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(probit_gp(y, dist, tau2, phi, iter));
    return rcpp_result_gen;
END_RCPP
}
// dmvnorm
double dmvnorm(arma::vec x, arma::vec mean, arma::mat sigma);
RcppExport SEXP _spmirt_dmvnorm(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnorm(x, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}
// testing
Rcpp::List testing(arma::mat X, arma::vec y);
RcppExport SEXP _spmirt_testing(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(testing(X, y));
    return rcpp_result_gen;
END_RCPP
}
// ifa_gibbs
Rcpp::List ifa_gibbs(Rcpp::NumericVector y, int n, int q, int N, int m);
RcppExport SEXP _spmirt_ifa_gibbs(SEXP ySEXP, SEXP nSEXP, SEXP qSEXP, SEXP NSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(ifa_gibbs(y, n, q, N, m));
    return rcpp_result_gen;
END_RCPP
}
// ifa_gibbs_nonide
Rcpp::List ifa_gibbs_nonide(Rcpp::NumericVector y, int n, int q, int N, int m);
RcppExport SEXP _spmirt_ifa_gibbs_nonide(SEXP ySEXP, SEXP nSEXP, SEXP qSEXP, SEXP NSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(ifa_gibbs_nonide(y, n, q, N, m));
    return rcpp_result_gen;
END_RCPP
}
// spifa_gibbs
Rcpp::List spifa_gibbs(Rcpp::NumericVector y, int n, int q, int N, int m);
RcppExport SEXP _spmirt_spifa_gibbs(SEXP ySEXP, SEXP nSEXP, SEXP qSEXP, SEXP NSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(spifa_gibbs(y, n, q, N, m));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _spmirt_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpptn_hello_world
List rcpptn_hello_world();
RcppExport SEXP _spmirt_rcpptn_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpptn_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spmirt_vec2mat", (DL_FUNC) &_spmirt_vec2mat, 3},
    {"_spmirt_vec2matt", (DL_FUNC) &_spmirt_vec2matt, 3},
    {"_spmirt_vecsub", (DL_FUNC) &_spmirt_vecsub, 3},
    {"_spmirt_vecsub1", (DL_FUNC) &_spmirt_vecsub1, 2},
    {"_spmirt_matsub1", (DL_FUNC) &_spmirt_matsub1, 3},
    {"_spmirt_subset_cpp", (DL_FUNC) &_spmirt_subset_cpp, 2},
    {"_spmirt_probit_gp", (DL_FUNC) &_spmirt_probit_gp, 5},
    {"_spmirt_dmvnorm", (DL_FUNC) &_spmirt_dmvnorm, 3},
    {"_spmirt_testing", (DL_FUNC) &_spmirt_testing, 2},
    {"_spmirt_ifa_gibbs", (DL_FUNC) &_spmirt_ifa_gibbs, 5},
    {"_spmirt_ifa_gibbs_nonide", (DL_FUNC) &_spmirt_ifa_gibbs_nonide, 5},
    {"_spmirt_spifa_gibbs", (DL_FUNC) &_spmirt_spifa_gibbs, 5},
    {"_spmirt_rcpp_hello_world", (DL_FUNC) &_spmirt_rcpp_hello_world, 0},
    {"_spmirt_rcpptn_hello_world", (DL_FUNC) &_spmirt_rcpptn_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_spmirt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
