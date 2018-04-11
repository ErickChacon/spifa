
library(MCMCpack)
library(purrr)

p <- 3
S <- diag(p)
n <- 10

riwish_tri <- function (v, S) {
  mat <- riwish(v, S)
  low_mat <- mat[lower.tri(mat, diag = TRUE)]
  return(as.numeric(low_mat))
}

panel.hist <- function(x, ...) {
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
   h <- hist(x, plot = FALSE)
   breaks <- h$breaks; nB <- length(breaks)
   y <- h$counts; y <- y/max(y)
   rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

cov_sim <- map(1:500, ~ riwish_tri(p + 8, S)) %>%
  reduce(rbind)
pairs(cov_sim, diag.panel = panel.hist)

summary(cov_sim)


coriwish_tri <- function (v, S) {
  mat <- riwish(v, S)
  mat_sd <- diag(diag(mat)^-0.5)
  mat <- mat_sd %*% mat %*% mat_sd
  # low_mat <- mat[lower.tri(mat, diag = TRUE)]
  low_mat <- mat[lower.tri(mat)]
  return(as.numeric(low_mat))
}

cor_sim <- map(1:500, ~ coriwish_tri(p + 1, S)) %>%
  reduce(rbind)
pairs(cor_sim, diag.panel = panel.hist)

coriwish_tri(3+1, S)

(bla <- riwish(3, S))
(var <- diag(diag(bla)^-0.5))
var %*% bla %*% var

mat <- matrix(1:3, nrow = 3, ncol = 3)
diag_mat <- diag(mat)

mat * diag_mat
t(mat) * diag_mat
# diag_mat * mat


p <- 3
A <- matrix(rnorm(p ^ 2), p)
(S <- crossprod(A))

X <- riwish(p+1, S)

sum(diag(S %*% solve(X)))
sum(S * solve(X))



library(dplyr)

test <- function (a, b) {
  a + b
}

test(3, 4)
library()

