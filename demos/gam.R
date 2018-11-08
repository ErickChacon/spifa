
setwd("../src")
getwd()

p <- 50
n <- 10
A <- matrix(rnorm(p ^ 2), p)
(S <- crossprod(A))
X <- MCMCpack::riwish(p+1, S)

sum(diag(S %*% solve(X)))
sum(S * solve(X))


Rcpp::sourceCpp("inv-wishart.cpp")
dinvwish(p+1, X, S)
MCMCpack::diwish(X, p+1, S)


require(rbenchmark)
benchmark(
  dinvwish(p + 1, X, S),
  MCMCpack::diwish(X, p + 1, S),
  columns = c('test', 'replications', 'relative', 'elapsed'),
  order = 'elapsed',
  replications = 3000
  )

