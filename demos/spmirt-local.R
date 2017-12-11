setwd("..")
getwd()

Rcpp::compileAttributes(getwd())
roxygen2::roxygenize(getwd(), roclets = c("collate", "namespace", "rd"))
path_repos <- system("find ~ -name Repositories 2>/dev/null", intern = T)
devtools::install_local(file.path(path_repos, "spmirt"))

# devtools::check()
# devtools::use_travis()


# setwd("demos")

library(spmirt)
library(Matrix)

rcpptn_hello_world()

set.seed(1)
(bla <- ifa_gibbs(rep(0, 500), 5, 5, 2))
hist(bla$z, 30)
# hist(bla$theta, 30)

bla <- ifa_gibbs(sample(0:1, 500, TRUE), 100, 5, 1)
# hist(bla$theta)

plot(bla$theta[, 1], type = "l")
hist(bla$theta[, 1], 100)
# matplot(bla$c[], type = "l")

bla$theta
bla$a
bla$c
bla$y
bla$z
bla$mu_z

plot(bla$z, bla$mu_z)

test(1:5)

X_theta(rep(1:3,5), 3, 5)

ls("package:spmirt")
# function 'RcppTN_rtn1' not provided by package 'RcppTN'

library(truncnorm)
library(RcppTN)

rtruncnorm(5, a=-Inf, b=0, mean = 0, sd = 1)
rcpptn_hello_world()
‘fill’ has not been declared


