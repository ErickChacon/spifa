
#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
Rcpp::StringVector name_samples_vec(int n_elem, std::string name) {
  Rcpp::StringVector param_names(n_elem);
  for (int i = 0; i < param_names.size(); ++i) {
    param_names(i) = name + "[" + std::to_string(i+1) + "]";
  }
  return param_names;
}

//' @export
// [[Rcpp::export]]
Rcpp::StringVector name_samples_mat(int nrow, int ncol, std::string name) {

  // determine the length of the output
  int n_elem = double(nrow * ncol);

  // fill the vector string
  int counter = 0;
  Rcpp::StringVector param_names(n_elem);
  for (int j = 0; j < ncol ; ++j) {
    for (int i = 0; i < nrow; ++i) {
      param_names(counter) = name + "[" + std::to_string(i+1) + "," +
        std::to_string(j+1) + "]";
      counter++;
    }
  }

  return param_names;
}

//' @export
// [[Rcpp::export]]
Rcpp::StringVector name_samples_lower(int nrow, int ncol, std::string name,
    bool diag = true) {

  // determine the length of the output
  int n_elem;
  if (nrow < ncol) {
    n_elem = round(nrow * (nrow + 1)) / 2;
  } else {
    n_elem = round(nrow * ncol - (ncol-1) * ncol / 2);
  }

  if (!diag) {
    n_elem = round(n_elem - std::min(nrow, ncol));
  }

  // fill the vector string
  int counter = 0;
  Rcpp::StringVector param_names(n_elem);
  for (int j = 0; j < ncol ; ++j) {
    for (int i = j + !diag; i < nrow; ++i) {
      param_names(counter) = name + "[" + std::to_string(i+1) + "," +
        std::to_string(j+1) + "]";
      counter++;
    }
  }

  return param_names;
}

