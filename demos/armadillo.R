library(spmirt)
x <- matrix(1:10, 5, 2)
x <- matrix(0, 7, 5)
(out <- subset_cpp(x, 1:10))
out$vec_scal
out$X

(out <- subset_cpp(x, 4:8))
out$vec_scal

(out <- subset_cpp(x, 8:4))
out$vec_scal

(vecsub(1:10, 1, 3))
# crossprod()

library(spmirt)
rcpptn_hello_world()

library(mvtnorm)

sigma <- matrix(c(1,0,0,0), ncol=2)
x <- rmvnorm(n=5, mean=c(1,2), sigma=sigma)
colMeans(x)
var(x)

dmvnorm(x=c(1, 2), mean = c(1,2), sigma = sigma)
dnorm(1, 1, 10^(-309))
dmvnorm(x=c(0,0), mean=c(1,1))

library(spmirt)
x <- matrix(0, 7, 5)
testing(x, 1:10)
