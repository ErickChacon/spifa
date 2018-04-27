// g++ armadillo-test.cpp -larmadillo -lopenblas
#include <iostream>
#include <armadillo>

int main(){
  arma::arma_version ver;
  std::cout << "ARMA version: "<< ver.as_string() << std::endl;
  arma::mat A = arma::randn<arma::mat>(5,5);
  for (int i = 0; i < 5; ++i) {
    A(i,i) = i+1;
  }
  // A.diag() = arma::span(1,5);
  // arma::mat U = arma::chol(A);
  // std::cout << arma::eye(5, 2) << std::endl;
  arma::rowvec indices(5, arma::fill::ones);
  indices(1) = 0;
  indices(3) = 0;
  std::cout << find(indices == 0) << std::endl;
  A.rows(find(indices == 0)) *= 0;
  A.cols(find(indices == 0)) *= 0;
  A.diag() += 1;
  std::cout << A << std::endl;
  std::cout << indices(2) << std::endl;
  std::cout << arma::join_cols(arma::ones(3), arma::zeros(4)) << std::endl;
  std::cout << indices.subvec(0, 4) << std::endl;

}
