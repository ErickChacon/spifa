#include <RcppArmadillo.h>

// [[Rcpp::export]]
std::string stringR(std::string text)
{
  if (text == "hello") {
    return "yes";
  } else {
    return "no";
  }
}

// [[Rcpp::export]]
arma::mat vecR(arma::mat input)
{
return input;
}
