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


library(mvtnorm)
library(microbenchmark)
sigma <- matrix(c(2,1,1,2), ncol=2)
z <- c(1.4, 5)
microbenchmark(
               mvtnorm::dmvnorm(z, mean=c(1,2), sigma=sigma, log = TRUE),
               spmirt::dmvnorm(z, mean=c(1,2), sigma=sigma)
)

colMeans(x)
var(x)

n <- 1000
f <- list(
  prob ~ I(1) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.02)),
  size ~ I(1)
  )
data <- sim_model(formula = f,
                  link_inv = list(pnorm, identity),
                  generator = rbinom,
                  n = n,
                  seed = 1)
data <- data %>% dplyr::arrange(s1)
dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
Sigma <- 2 * exp(-dist/0.02)
z <- rnorm(n)
microbenchmark(
  mvtnorm::dmvnorm(z, mean=rep(0, n), sigma=Sigma, log = TRUE),
  spmirt::dmvnorm(z, mean=rep(0, n), sigma=Sigma),
  times = 10
)

# mvtnorm::dmvnorm(z, mean=rep(0, n), sigma=Sigma, log = TRUE)
# spmirt::dmvnorm(z, mean=rep(0, n), sigma=Sigma)

psych::logistic(2)
psych::logit(0.5)


hist(psych::logit(runif(10000)), 100)
plop <- runif(10000)
hist(plop / (1-plop), 100)

