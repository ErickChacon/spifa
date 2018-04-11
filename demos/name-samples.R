

getwd()
Rcpp::sourceCpp("../src/names.cpp")

name_samples_vec(5, "beta")

name_samples_lower(6, 3, "corr")
name_samples_lower(3, 6, "corr")

name_samples_lower(6, 3, "corr", FALSE)
name_samples_lower(3, 6, "corr", FALSE)

name_samples_mat(3,2, "sigma")



