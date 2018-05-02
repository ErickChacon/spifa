// g++ reference-test.cpp -larmadillo -lopenblas
#include <iostream>
#include <armadillo>

void change(const arma::mat& A)
{
  // A.diag().fill(3);
  // A.diag() += 3;
  arma::mat B = A + 1;
  std::cout << B << std::endl;
}

int main(){
  const arma::mat A = arma::randn<arma::mat>(5,5);
  // for (int i = 0; i < 5; ++i) {
  //   A(i,i) = i+1;
  // }
  std::cout << "Testing diagonal row operator:" << std::endl;
  std::cout << A << std::endl;
  change(A);
  std::cout << A << std::endl;

}

