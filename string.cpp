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
arma::mat vecR(arma::mat input)
{
return input;
}


// [[Rcpp::export]]
double nullR(double plop)
{
  return plop;
}
