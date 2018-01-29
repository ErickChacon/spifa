library(spmirt)
library(ggplot2)
library(datasim)

n <- 300

f <- list(
  prob ~ I(0) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.02),
       sigma2 = 5),
  size ~ I(1)
  )
data <- sim_model(formula = f,
                  link_inv = list(pnorm, identity),
                  generator = rbinom,
                  n = n)
data <- dplyr::rename(data, gp = gp.list.prob)

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = gp, size = gp))

vg <- gstat::variogram(gp ~ 1, ~ s1 + s2, data, cutoff = 0.4, width = 0.005)
ggplot(vg, aes(dist, gamma, weight = np)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.4))

vg <- gstat::variogram(response ~ 1, ~ s1 + s2, data, cutoff = 0.4, width = 0.001)
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.4))


iter <- 1000
dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
out <- probit_gp(data$response, dist, c(psych::logit(0.5), log(0.02)), iter)
plot(out$param)

data$z <- out$z[iter,]
ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = gp, size = gp))

vg <- gstat::variogram(z ~ 1, ~ s1 + s2, data, cutoff = 0.4, width = 0.005)
ggplot(vg, aes(dist, gamma, weight = np)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.4))

hist(out$z[iter,])
# hist(out$z[2,])
# hist(out$z[3,])
# hist(out$z[4,])
# hist(out$z[5,])
