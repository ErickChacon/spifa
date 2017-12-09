setwd("..")
getwd()

Rcpp::compileAttributes(getwd())
roxygen2::roxygenize(getwd(), roclets = c("collate", "namespace", "rd"))
path_repos <- system("find ~ -name Repositories 2>/dev/null", intern = T)
devtools::install_local(file.path(path_repos, "spmirt"))

devtools::check()
devtools::use_travis()

library(spmirt)


