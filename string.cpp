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
std::vector<std::string> string_in(std::vector<std::string> v)
{
  if (std::find(v.begin(), v.end(), "c") != v.end())
  {
    // Element in vector.
    Rcpp::Rcout << "it contains" << std::endl;
  }

  // Rcpp::StringVector myclass(2);
  // myclass = {"a", "b"};
  // Rcpp::Rcout << myclass << std::endl;
  // myclass(0) = "spmirt.list";
  // myclass(1) = "list";
  // output.attr("class") = myclass;


  // std::vector<std::string> text;
  // text = {"one", "two", "three"};
  // if (text == "hello") {
  //   return "yes";
  // } else {
  //   return "no";
  // }
  return v;
}

// [[Rcpp::export]]
Rcpp::StringVector string_subset(Rcpp::StringVector v)
{

  Rcpp::StringVector out;
  arma::mat T = arma::eye(4, 3);
  T(1,0) = 1;
  arma::uvec index = find(T == 1);
  Rcpp::Rcout << index << std::endl;
  // Rcpp::IntegerVector idx = Rcpp::IntegerVector::create(0, 1, 2);
  // Rcpp::IntegerVector idx = Rcpp::wrap(index);
  Rcpp::IntegerVector idx = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(index));
  out = v[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(index))];

  return out;
}




// [[Rcpp::export]]
arma::mat vecR(arma::mat input)
{
return input;
}


// [[Rcpp::export]]
double nullR(double plop)
{
  return plop;
}
