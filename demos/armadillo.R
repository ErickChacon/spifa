library(spmirt)
x <- matrix(1:10, 5, 2)
(out <- subset_cpp(x, 1:10))
# lapply(out, class)
