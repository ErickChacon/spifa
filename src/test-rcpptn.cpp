
#include <Rcpp.h>
#include <RcppTN.h>
// [[Rcpp::depends(RcppTN)]]

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
List rcpptn_hello_world() {
  double a = RcppTN::rtn1(0.0, 1.0, R_NegInf, 3.7) ;
  double a2 = RcppTN::rtn1(0.0, 1.0, R_NegInf, 3.7) ;
  double b = RcppTN::etn1(0.0, 1.0, 3.5, 3.7) ;
  double c = RcppTN::vtn1(0.0, 1.0, 3.5, 3.7) ;
  double d = RcppTN::dtn1(3.6, 0.0, 1.0, 3.5, 3.7) ;
  double e = RcppTN::enttn1(0.0, 1.0, 3.5, 3.7) ;
  double cst = RcppTN::rtn1(1.0, 0.0, 0, 3.7) ;
  NumericVector y = NumericVector::create(a, a2, b, c, d, e) ;
  List z = List::create( y,
      cst) ;
  return(z) ;
}
