
library(tidyverse)
getwd()
# Rcpp::sourceCpp("../src/Power-driver.cpp")
# classtest(4, 3)
Rcpp::sourceCpp("../src/ifa-main.cpp")
# Rcpp::sourceCpp("../src/ifa-driver.cpp")
(samples <- spmirt(response = sample(0:1, 10, TRUE), nobs = 5, nitem = 2, nfactors = 1,
                   L_rest = matrix(1, 2, 1), niter = 1000))

