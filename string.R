Rcpp::sourceCpp("string.cpp")
stringR("hellod")
vecR(matrix(NA))
string_in(letters[1:12])
string_subset(letters[1:12])

matrix(letters[1:12], 4, 3)

x <- 1:10
x[sample(1:length(x), 4)] <- NA
na_fill(c(0,1,NA))

R <- matrix(1, 3, 3)
# diag(R) = 1
DRD(R, 1:3)
diag(1:3) %*% R %*% diag(1:3)

R <- matrix(rnorm(9), 3)
R <- R %*% t(R)
(bla <- t(chol(R)))
solve(bla)
