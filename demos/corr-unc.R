
# library(spmirt)

getwd()
Rcpp::sourceCpp("../src/correlation.cpp")


k <- 3
# N <- k * (k-1)/2
N <- k * (k+1)/2
x <- 1:N
(bla <- vec2trimatl(x, 3, TRUE))
(x2 <- trimatl2vec(bla, TRUE))

N <- k * (k-1)/2
x <- 1:N
(bla <- vec2trimatl(x, 3, FALSE))
(x2 <- trimatl2vec(bla, FALSE))

N <- k * (k-1)/2
x <- rnorm(N)
(bla <- vec2chol_corr(x, 3))
(bla2 <- chol_corr2vec(bla))

tcrossprod(bla)
# t(chol(mycorr))


