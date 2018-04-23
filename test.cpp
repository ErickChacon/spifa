#include <iostream>
#include <R.h>
#include <Rinternals.h>

double randomise() {
  std::cout << "hello\n";
  // std::cout << R::unif(0,1) << "\n";
  double out = rnorm(0.0,1.0);
  return out;
}

double function_name(double x) {
  x += 1;
  return x;
}

