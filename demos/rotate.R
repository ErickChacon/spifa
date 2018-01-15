rotate <- function (theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
}
rotate(pi/2)

aliase <- function(signs = c(1, 1), switch = TRUE) {
  out <- diag(signs)
  if (switch) {
    out = out[, 2:1]
  }
  return(out)
}

# trans <- matrix(1, 2, 2)
(trans <- rotate(pi/2))
(trans <- rotate(pi/3))
t(trans) %*% diag(c(1, 1.3)) %*% trans
t(trans) %*% trans

rotate(0)
rotate(pi/2)
rotate(pi)
rotate(3*pi/2)

point = matrix(c(0.5, 1.5))
plot(t(point), pch = 19, xlim = c(-2,2), ylim = c(-2, 2))
points(t(rotate(pi/2) %*% point), pch = 19, col = 2)
points(t(rotate(pi) %*% point), pch = 19, col = 3)
points(t(rotate(3*pi/2) %*% point), pch = 19, col = 4)
points(t(rotate(3*pi/2) %*% point), pch = 19, col = 4)
points(t(aliase() %*% point), pch = 17, col = 6)
points(t(aliase(c(-1, -1)) %*% point), pch = 17, col = 6)
points(t(aliase(c(-1, 1), FALSE) %*% point), pch = 17, col = 5)
points(t(aliase(c(1, -1), FALSE) %*% point), pch = 17, col = 5)
abline(v = c(3/4, 4/5))

t(rotate(pi/4)) %*% rotate(pi/4)
t(rotate(pi/4)) %*% matrix(c(-1, 0, 0, 1), 2) %*% rotate(pi/4)

rotate(pi/4) %*% rbind(1, 1)



