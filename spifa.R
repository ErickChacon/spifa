
Rcpp::compileAttributes(".", verbose = TRUE)                                         # cpp to R
devtools::load_all(".")                                              # compile src code before document
roxygen2::roxygenize(".", roclets = c("collate", "namespace", "rd")) # document

devtools::check_built(".")



devtools::load_all(".")

