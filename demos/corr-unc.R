
# library(spmirt)

getwd()
Rcpp::sourceCpp("../src/correlation.cpp")


k <- 3
# N <- k * (k-1)/2
N <- k * (k+1)/2
x <- 1:N
(bla <- vec2trimatl(x, 3, TRUE))
class(bla)

N <- k * (k-1)/2
x <- 1:N
(bla <- vec2trimatl(x, 3, FALSE))
class(bla)

N <- k * (k-1)/2
x <- rnorm(N)
(bla <- vec2corr(x, 3))
class(bla)

tcrossprod(bla)
# t(chol(mycorr))


