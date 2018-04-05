
# library(spmirt)

getwd()
Rcpp::sourceCpp("../src/correlation.cpp")


k <- 5
N <- k * (k+1)/2
x <- 1:N
(trimatl_x <- vec2trimatl(x, k, TRUE))
x2 <- trimatl2vec(trimatl_x, TRUE)
(as.numeric(x2))
(x)

N <- k * (k-1)/2
x <- 1:N
(trimatl_x <- vec2trimatl(x, k, FALSE))
x2 <- trimatl2vec(trimatl_x, FALSE)
(as.numeric(x2))
(x)

N <- k * (k-1)/2
x <- rnorm(N)
(chol_corr <- vec2chol_corr(x, k))
x2 <- chol_corr2vec(chol_corr)
(as.numeric(x2))
(x)

Corr <- tcrossprod(chol_corr)


getwd()
Rcpp::sourceCpp("../src/likelihood.cpp")

k <- 2
N <- k * (k-1)/2
x <- rnorm(N)
chol_corr <- vec2chol_corr(x, k)
(Corr <- tcrossprod(chol_corr))
dlkj_corr(Corr, 2)

cor_seq <- seq(-0.99, 0.99, length.out = 100)
dcor_seq <- purrr::map_dbl(cor_seq, ~ dlkj_corr(matrix(c(1, ., ., 1), 2), 10.5))
plot(cor_seq, dcor_seq)

cor_chol_seq <- seq(-1, 1, length.out = 100)
dcor_chol_seq <- purrr::map_dbl(cor_chol_seq,
                           ~ dlkj_corr_chol(vec2chol_corr(., 2), 1))
plot(cor_chol_seq, dcor_chol_seq)

tcrossprod(vec2chol_corr(-0.5, 2))

vec2chol_corr(rnorm(1), 2)


getwd()
Rcpp::sourceCpp("../src/likelihood.cpp")

k <- 3
N <- k * (k-1)/2
x <- rnorm(N)
(Corr <- tcrossprod(vec2chol_corr(x, k)))
dlkj_corr(Corr, 2)

seq <- seq(-0.99, 0.99, length.out = 40)
cor_seq <- expand.grid(z1 = seq, z2 = seq)
dcor_seq <- purrr::map_dbl(cor_seq,
               ~ dlkj_corr(tcrossprod(vec2chol_corr(c(.$z1, $z2), k)), 5))

plot(cor_seq, dcor_seq)

cor_chol_seq <- seq(-1, 1, length.out = 100)
dcor_chol_seq <- purrr::map_dbl(cor_chol_seq,
                           ~ dlkj_corr_chol(vec2chol_corr(., 2), 1))
plot(cor_chol_seq, dcor_chol_seq)
