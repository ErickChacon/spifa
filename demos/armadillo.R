library(spmirt)
x <- matrix(1:10, 5, 2)
(out <- subset_cpp(x, 1:10))
out$vec_scal

(out <- subset_cpp(x, 4:8))
out$vec_scal

(out <- subset_cpp(x, 8:4))
out$vec_scal

(vecsub(1:10, 1, 3))
# crossprod()

