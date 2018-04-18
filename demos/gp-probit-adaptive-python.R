library(spmirt)
library(ggplot2)
library(datasim)

n <- 300

f <- list(
  prob ~ I(0) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.04),
       sigma2 = 1),
  size ~ I(1)
  )

data <- sim_model(formula = f,
                  link_inv = list(pnorm, identity),
                  generator = rbinom,
                  n = n,
                  seed = 2
                  )
data <- dplyr::rename(data, gp = gp.list.prob)

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = gp, size = gp)) +
  scale_colour_distiller(palette = "RdYlBu")

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = factor(response)), size = 2)

data$gp2 = data$gp + rnorm(300)


vg <- gstat::variogram(gp2 ~ 1, ~ s1 + s2, data, cutoff = 0.7, width = 0.005)
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

noise <- rnorm(n)


# library(reticulate)
# spy <- import_from_path("spmirtpy.models.gpreg",
#                         path = "/home/chaconmo/Documents/Repositories")
reticulate::source_python("~/Documents/Repositories/spmirtpy/models/gpreg.py")

iter <- 5000L
dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
# sigma_prop <- matrix(c(0.1, -0.02, -0.02, 0.05), 2) * 2.38 ^ 2 / 2
# sigma_prop <- matrix(c(0.2, -0.02, -0.02, 0.05), 2) * 2.38 ^ 2 / 2
# sigma_prop <- matrix(c(0.2, -0.02, -0.02, 0.05), 2) * 2.38 ^ 2 / 2
sigma_prop0 <- matrix(c(0.034, 0.011, 0.011, 0.067), 2) * 2.38 ^ 2 / 2
sigma_prop <- matrix(c(0.034, 0.011, 0.011, 0.067), 2) * 1.65 ^ 2 / (2 ^ (1/3))
# sigma_prop <- matrix(c(0.1, 0.04, 0.04, 0.13), 2) * 1.65 ^ 2 / (2 ^ (1/3))
# sigma_prop <- matrix(c(0.27, -0.07, -0.07, 0.21), 2) * 2.38 ^ 2 / 2
# system.time(
#   (plop <- gpreg(cbind(data$gp + rnorm(n)) , dist, cbind(c(log(1.0), log(0.04))), iter, sigma_prop))
#   )
# plot(plop[,1])
# plot(plop[,2])
# # plop
# nrow(unique(plop)) / iter
# #   print(plop)
# # plop[1:5, 1:5]
# # test()
# plop <- geo_normal(cbind(data$gp + noise) , dist, cbind(c(log(0.95), log(0.07))))
plop <- geo_normal(cbind(data$gp + noise) , dist, cbind(c(log(0.95), log(0.03))))
plop$map()


# system.time(
# plop1 <- plop$sample(sigma_prop0, iter, "rwmh")
# )
# system.time(
# plop2 <- plop$sample(sigma_prop, iter, "mala", 1)
# )
#
# opar <- par(mfrow = c(1, 2))
# plot(density(plop1[,1]))
# lines(density(plop2[,1]), col = 2)
# plot(density(plop1[,2]))
# lines(density(plop2[,2]), col = 2)
# par(opar)


# plot(plop1[,1])
# points(plop2[,1], col = 2)
#
# plot(plop1[,2], type = "l")
# lines(plop2[,2], type = "l", col = 2)
# abline(h = 0.04, col = 2)
# plop
nrow(unique(plop1)) / iter
nrow(unique(plop2)) / iter
summary(plop2[,2])

coda::effectiveSize(plop1)
coda::effectiveSize(plop2)


# initialize(plop$callpost(), plop$feed_dict)

# plop$posterior(cbind(c(log(1.0), log(0.04))), 1)

# dmvnorm(c(0,0), cbind(c(0,0)), matrix(c(1,0.5,0.5,2)))
# noise = rnorm(n)
# test(cbind(c(log(1.0), log(0.04))), dist, cbind(data$gp + noise))



Rcpp::sourceCpp("../src/gp-gibss-adap.cpp")

# 0.6
# 50000 -> 17 minutes
iter <- 10000
dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
# out <- probit_gp(data$response, dist, c(psych::logit(0.5), log(0.02)), iter)
# sigma_prop <- matrix(c(0.1, 0.05, 0.05, 0.1), 2) / 10
# sigma_prop <- matrix(c(0.1, -0.02, -0.02, 0.05), 2) * 2.38^2/2
# sigma_prop <- matrix(c(0.1, -0.02, -0.02, 0.05), 2) * 2.38 ^ 2 / 2
# system.time(
#   out2 <- probit_gp_adap(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )
system.time(
  out <- probit_gp_chol2(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
)
# system.time(
#   out0 <- probit_gp(data$response, dist, c(log(1), log(0.05)), iter, sigma_prop)
# )

# plot(out$param[, 2])
# abline(h = 0.03, col = 2)

beta <- 0.05
x <- 1
y <- -1

mixture <- function (x, y) {
  (1-beta) * dnorm(y, x, 2.38 * 0.5) + beta * dnorm(y, x, 0.1)
}

mixture(x, y)
mixture(y, x)


plot(out$param, type = "b")
points(1, 0.04, col = 2, pch = 19)

plot(out$param[(iter/2):iter,], type = "b")
points(1, 0.04, col = 2, pch = 19)

plot(out$param[seq(iter/2, iter, 10), ], type = "b")
points(1, 0.04, col = 2, pch = 19)

df <- setNames(as.data.frame(out$param[seq(iter/2, iter, 15), ]), c("sigma2", "phi"))
ggplot(data = df, aes(sigma2, phi)) +
  # stat_density2d() +
  stat_density2d(aes(fill=..level..,alpha=..level..),
                 geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  # geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  guides(alpha="none") +
  geom_point(alpha = 0.5)

df2 <- setNames(as.data.frame(log(out$param[seq(iter/2, iter, 15), ])), c("sigma2", "phi"))
ggplot(data = df2, aes(sigma2, phi)) +
  # stat_density2d() +
  stat_density2d(aes(fill=..level..,alpha=..level..),
                 geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  # geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  guides(alpha="none") +
  geom_point(alpha = 0.5)



nrow(unique(out$param))/(iter)

nrow(unique(out$param[(iter/2 + 1):iter, ]))/(iter/2)

plot(out$param[, 1], type = "b")
abline(h = 1, col = 2)

plot(out$param[seq(iter/ 2, iter, 10), 1], type = "b")
abline(h = 1, col = 2)

plot(out$param[, 2], type = "b")
abline(h = 0.04, col = 2)
summary(out$param[, 2])

plot(out$param[seq(iter/2, iter, 10), 2], type = "b")
abline(h = 0.04, col = 2)
summary(out$param[seq(iter/2, iter, 10), 2])

# var(out$param)
var(out$param[(iter/2):iter,])
var(log(out$param[(iter/2):iter,]))

# varfun1 <- function (i) {
#   var(log(out$param[i + 1:1000,1]))
# }
# varfun2 <- function (i) {
#   var(log(out$param[i + 1:1000,2]))
# }
# varfun3 <- function (i) {
#   cov(log(out$param[i + 1:1000, 1]), log(out$param[i + 1:1000, 2]))
# }
varfun1 <- function (i) {
  var(log(out$param[1:i,1]))
}
varfun2 <- function (i) {
  var(log(out$param[1:i,2]))
}
varfun3 <- function (i) {
  cov(log(out$param[1:i, 1]), log(out$param[1:i, 2]))
}
index <- seq(3, iter, 10)
plot(index, sapply(index, varfun1), type = "b")

plot(index, sapply(index, varfun2), type = "b")

plot(index, sapply(index, varfun3), type = "b")


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


# # phi
# hist(exp(rnorm(10000, log(0.03), 0.4)), 100)
# hist(rnorm(10000, log(0.03), 0.4), 100)
#
# # sigma^2
# hist(exp(rnorm(10000, log(1), 0.4)), 100)
# hist(rnorm(10000, log(1), 0.4), 100)

# variance_lognormal <- function (sigma2, mu) (exp(sigma2) - 1) * exp(2 * mu + sigma2)
# sigma2 <- seq(0.001, 5, length.out = 100)
# plot(sigma2, variance_lognormal(sigma2, 0.03))
#

# phi
hist(exp(rnorm(50000, log(0.03), 0.4)), 200, freq = FALSE )
lines(density(out$param[(iter/2):iter, 2]), col = 2)
quantile(out$param[(iter/2):iter, 2], c(0.025, 0.975))

# sigma^2
hist(exp(rnorm(50000, log(1), 0.4)), 200, freq = FALSE, xlim = c(0, 5))
# hist(exp(rnorm(5000, log(1), 0.3)), 100, freq = FALSE)
lines(density(out$param[(iter/2):iter, 1]), col = 2)
quantile(out$param[(iter/2):iter, 1], c(0.025, 0.975))

# acf(out$param[(iter/2):iter, 1])
# acf(out$param[(iter/2):iter, 2])

acf(out$param[seq((iter/2), iter, 10), 1])
# acf(out$param[seq((iter/2), iter, 10), 2])

index <- seq(0, iter - 1000, 50)
plot(index, sapply(index, varfun1), type = "b")

plot(index, sapply(index, varfun2), type = "b")

plot(index, sapply(index, varfun3), type = "b")


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


# # phi
# hist(exp(rnorm(10000, log(0.03), 0.4)), 100)
# hist(rnorm(10000, log(0.03), 0.4), 100)
#
# # sigma^2
# hist(exp(rnorm(10000, log(1), 0.4)), 100)
# hist(rnorm(10000, log(1), 0.4), 100)

# variance_lognormal <- function (sigma2, mu) (exp(sigma2) - 1) * exp(2 * mu + sigma2)
# sigma2 <- seq(0.001, 5, length.out = 100)
# plot(sigma2, variance_lognormal(sigma2, 0.03))
#

# phi
hist(exp(rnorm(50000, log(0.03), 0.4)), 200, freq = FALSE )
lines(density(out$param[(iter/2):iter, 2]), col = 2)
quantile(out$param[(iter/2):iter, 2], c(0.025, 0.975))

# sigma^2
hist(exp(rnorm(50000, log(1), 0.4)), 200, freq = FALSE, xlim = c(0, 5))
# hist(exp(rnorm(5000, log(1), 0.3)), 100, freq = FALSE)
lines(density(out$param[(iter/2):iter, 1]), col = 2)
quantile(out$param[(iter/2):iter, 1], c(0.025, 0.975))

acf(out$param[(iter/2):iter, 1])
acf(out$param[(iter/2):iter, 2])

acf(out$param[seq((iter/2), iter, 10), 1])
# acf(out$param[seq((iter/2), iter, 10), 2])

