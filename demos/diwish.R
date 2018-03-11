function (W, v, S) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop("S not square in diwish().\n")
    }
    if (!is.matrix(W)) 
        W <- matrix(W)
    if (nrow(W) != ncol(W)) {
        stop("W not square in diwish().\n")
    }
    if (nrow(S) != ncol(W)) {
        stop("W and X of different dimensionality in diwish().\n")
    }
    if (v < nrow(S)) {
        stop("v is less than the dimension of S in  diwish().\n")
    }
    p <- nrow(S)
    gammapart <- sum(lgamma((v + 1 - 1:p)/2))
    ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 
        1) * log(pi)
    cholS <- chol(S)
    cholW <- chol(W)
    halflogdetS <- sum(log(diag(cholS)))
    halflogdetW <- sum(log(diag(cholW)))
    invW <- chol2inv(cholW)
    exptrace <- sum(S * invW)
    lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * 
        exptrace
    lpdf <- lnum - ldenom
    return(exp(lpdf))
}


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

