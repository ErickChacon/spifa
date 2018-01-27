library(spmirt)
library(ggplot2)
library(datasim)

n <- 1000

f <- list(
  prob ~ I(1) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.02)),
  size ~ I(1)
  )

data <- sim_model(formula = f,
                  link_inv = list(pnorm, identity),
                  generator = rbinom,
                  n = n,
                  seed = 1)
data <- dplyr::rename(data, gp = gp.list.prob)

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = gp, size = gp))

vg <- gstat::variogram(gp ~ 1, ~ s1 + s2, data, cutoff = 0.4, width = 0.005)
ggplot(vg, aes(dist, gamma, weight = np)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.4))

iter <- 20
dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
out <- probit_gp(data$response, dist, 1, 0.02, iter)

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
