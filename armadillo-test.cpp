#include <iostream>
#include <armadillo>

int main(){
  arma::arma_version ver;
  std::cout << "ARMA version: "<< ver.as_string() << std::endl;
}
