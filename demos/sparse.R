

y <- data_long$y
n <- length(unique(data_long$subject))
q <- length(unique(data_long$item))
x_c <- kronecker(Matrix(rep(1, n)), Diagonal(q))
x_theta <- 
image(x_c)

bla <- kronecker(Diagonal(2), Matrix(rep(1, 4), 2, 2))
print(bla)
image(bla)

     M1 <- Matrix(0+0:5, 2,3)
     (M <- kronecker(Diagonal(3), M1))

image(Diagonal(1, 4))

X_theta_c <- 
  sparseMatrix()
as()

plot(diag.spam(1, 9))

f <- gl(10, 2)

t(as(f, Class = "sparseMatrix"))

library(Matrix)

i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
A <- sparseMatrix(i, j, x = x)

print(A)

image(A)
