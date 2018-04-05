
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
Rcpp::sourceCpp("../src/pdf.cpp")

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
                           ~ dlkj_corr_chol(vec2chol_corr(., 2), 5))
plot(cor_chol_seq, dcor_chol_seq)



getwd()
Rcpp::sourceCpp("../src/pdf.cpp")

k <- 3
N <- k * (k-1)/2
x <- rnorm(N)
chol_corr <- vec2chol_corr(x, k)
(Corr <- tcrossprod(chol_corr))
dlkj_corr(Corr, 2)

cor_seq <- seq(-0.99, 0.99, length.out = 100)
dcor_seq <- purrr::map_dbl(cor_seq, ~ dlkj_corr(matrix(c(1, ., ., 1), 2), 1.5))
plot(cor_seq, dcor_seq)

cor_chol_seq <- seq(-1, 1, length.out = 100)
dcor_chol_seq <- purrr::map_dbl(cor_chol_seq,
                           ~ dlkj_corr_chol(vec2chol_corr(., 2), 1.5))
plot(cor_chol_seq, dcor_chol_seq)

cor_seq <- seq(-1.99, 1.99, length.out = 100)
dcor_seq <- purrr::map_dbl(cor_seq, ~ dlkj_corr_free(c(.), 2, 1.5))
plot(cor_seq, dcor_seq)

