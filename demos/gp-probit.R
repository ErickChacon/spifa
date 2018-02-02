library(spmirt)
library(ggplot2)
library(datasim)

n <- 1000

f <- list(
  prob ~ I(0) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.03),
       sigma2 = 1),
  size ~ I(1)
  )
data <- sim_model(formula = f,
                  link_inv = list(pnorm, identity),
                  generator = rbinom,
                  n = n
                  # ,
                  # seed = 2
                  )
data <- dplyr::rename(data, gp = gp.list.prob)

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = gp, size = gp)) +
  scale_colour_distiller(palette = "RdYlBu")

vg <- gstat::variogram(gp ~ 1, ~ s1 + s2, data, cutoff = 0.7, width = 0.005)
ggplot(vg, aes(dist, gamma, weight = np)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.7))

vg <- gstat::variogram(response ~ 1, ~ s1 + s2, data, cutoff = 0.7, width = 0.005)
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.7))


iter <- 10
dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
# out <- probit_gp(data$response, dist, c(psych::logit(0.5), log(0.02)), iter)
# sigma_prop <- matrix(c(0.1, 0.05, 0.05, 0.1), 2) / 10
sigma_prop <- matrix(c(0.1, 0, 0, 0.1), 2) / 10
# system.time(
#   out0 <- probit_gp(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )
# system.time(
#   out <- probit_gp_chol(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )
system.time(
  out <- probit_gp_chol2(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
)

# 34%
# DGEMV  performs one of the matrix-vector operations
#
#     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
#
#  where alpha and beta are scalars, x and y are vectors and A is an
#  m by n matrix.

# plot(out$param[, 2])
# abline(h = 0.03, col = 2)

plot(out$param, type = "b")
plot(out$param[, 1])
plot(out$param[, 2])
abline(h = 0.03, col = 2)
summary(out$param[, 2])
var(out$param)
# summary(out$param[, 2])

# data$z <- out$z[iter,]
# ggplot(data, aes(s1, s2)) +
#   geom_point(aes(col = gp, size = gp))
# vg <- gstat::variogram(z ~ 1, ~ s1 + s2, data, cutoff = 0.4, width = 0.005)
# ggplot(vg, aes(dist, gamma)) +
#   geom_point(aes(size = np)) +
#   geom_smooth() +
#   expand_limits(y = 0, x = 0) +
#   scale_x_continuous(limits = c(0, 0.4))
#
# hist(out$z[iter,])
# plot(out$z[iter,], data$response)
#
# plot(out$z[iter,], type = "b")
# var(out$z[iter,])
# # hist(out$z[2,])
# # hist(out$z[3,])
# # hist(out$z[4,])
# # hist(out$z[5,])

#
# # phi
# hist(exp(rnorm(1000, log(0.03), 0.4)))
#
# # sigma^2
# hist(exp(rnorm(1000, log(1), 0.4)))
