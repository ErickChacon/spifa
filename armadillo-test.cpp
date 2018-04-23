#include <iostream>
#include <armadillo>

int main(){
  arma::arma_version ver;
  std::cout << "ARMA version: "<< ver.as_string() << std::endl;
  arma::mat A = 0.5 * arma::eye<arma::mat>(5,5);
  arma::mat U = arma::chol(A);
  std::cout << U << std::endl;
}
