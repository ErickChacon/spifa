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
  std::cout << "Testing subsetting:" << std::endl;
  arma::vec low_thresh;
  low_thresh << 1 << 2 << arma::endr;
  std::cout << low_thresh << std::endl;
  std::cout << "Testing diagonal row operator:" << std::endl;
  std::cout << A << std::endl;
  A.diag() += 1 / A.row(0);
  std::cout << A << std::endl;
  A.diag() += A.col(0);
  std::cout << A << std::endl;
  arma::vec b;
  b << 1 << 2 << 3 << 4 << 5 << arma::endr;
  arma::vec b2 = b - 3;
  std::cout << b << std::endl;
  std::cout << b2 << std::endl;
  b %= b2;
  std::cout << b << std::endl;
  // low_thresh(0) = R_NegInf;
  // low_thresh(1) = 0;
  int n1 = 100;
  int n2 = 3;
  std::cout << 100/3 << std::endl;
  std::cout << 100 % 3 << std::endl;
  std::cout << 0 % 10 << std::endl;

  arma::mat X = arma::randn<arma::mat>(5,10);
  std::cout << arma::mean(arma::mean(X)) << std::endl;
  std::cout << arma::mean(arma::stddev(X)) << std::endl;
  std::cout << arma::mean(arma::stddev(X).t()) << std::endl;
  double ba = arma::as_scalar(arma::mean(arma::stddev(X).t()));

  double y = 6;
  y /= 2;
  std::cout << y << std::endl;

  std::cout << "Testing fill sparse" << std::endl;

  arma::mat T = arma::eye(4, 3);
  T(1,0) = 1;
  // T(2,1) = 1;
  std::cout << T << std::endl;
  std::cout << find(T != 0) << std::endl;
  arma::vec t = arma::regspace(1, 4);
  T(find(T != 0)) = t;
  std::cout << T << std::endl;

  std::cout << "Testing fill sparse by row" << std::endl;
  std::cout << T << std::endl;
  std::cout << find(T != 0) << std::endl;
  arma::uvec ix = find(T != 0);
  // std::cout << ix % 4.0 << std::endl;
  arma::umat im = arma::ind2sub(size(T), ix);
  std::cout << im << std::endl;
  std::cout << find(im.row(0) == 0) << std::endl;
  std::cout << t << std::endl;
  for (int i = 0; i < 4; ++i) {
    t(find(im.row(0) == i)) /= (i+1);
  }
  std::cout << t << std::endl;
  // arma::vec t = arma::regspace(1, 4);
  // T(find(T != 0)) = t;
  // std::cout << T << std::endl;
}
