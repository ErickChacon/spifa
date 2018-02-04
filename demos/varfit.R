
library(spmirt)
library(ggplot2)
library(datasim)

n <- 8000

f <- list(
  prob ~ I(0.5) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.02),
       sigma2 = 0.5),
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

ggplot(data, aes(s1, s2)) +
  geom_point(aes(col = factor(response)), size = 2)

vg0 <- gstat::variogram(gp ~ 1, ~ s1 + s2, data, cutoff = 0.7, width = 0.005)
ggplot(vg0, aes(dist, gamma, weight = np)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.7))

get_sill <- function (beta0, sigma2) {
  prop <- mean(pnorm(rnorm(10^5, beta0, sd = sigma2^0.5)))
  return(prop * (1 - prop))
}

get_nugget <- function (beta0, sigma2) {
  sill <- get_sill(beta0, sigma2)
  nugget <- sill - var(pnorm(rnorm(10^5, beta0, sd = sigma2^0.5)))
  return(nugget)
}

get_sigma2 <- function (beta0, sigma2) {
  sigma2 <- var(pnorm(rnorm(10^5, beta0, sd = sigma2^0.5)))
  return(sigma2)
}

vg <- gstat::variogram(response ~ 1, ~ s1 + s2, data, cutoff = 0.5, width = 0.005)
nugget <- get_nugget(0.5, 1)
sigma2 <- get_sigma2(0.5, 1)
phi_aux <- 0.02
vg$theor <- nugget + sigma2 * (1 - exp(-vg$dist / phi_aux))
inter <- 0.43
sig <- 0.048
phi_aux <- 0.047
vg$theor1 <- get_nugget(inter, sig) + get_sigma2(inter, sig) * (1 - exp(-vg$dist / phi_aux))
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  geom_line(aes(y = theor), col = 2) +
  geom_line(aes(y = theor1), col = 3) +
  scale_x_continuous(limits = c(0, 0.5)) +
  expand_limits(y = 0, x = 0)
library(gstat)
fit.variogram(vg, vgm("Exp"))

cost <- function (param, vg) {
  tau2 <- exp(param[1])
  sigma2 <- exp(param[2])
  phi <- exp(param[3])
  theor <- tau2 + sigma2 * (1 - exp(-vg$dist/phi))
  sum(vg$np * (theor - vg$gamma) ^ 2)
}
get_sill <- function (beta0, sigma2) {
  prop <- mean(pnorm(rnorm(10^5, beta0, sd = sigma2^0.5)))
  return(prop * (1 - prop))
}
get_nugget <- function (beta0, sigma2) {
  sill <- get_sill(beta0, sigma2)
  nugget <- sill - var(pnorm(rnorm(10^5, beta0, sd = sigma2^0.5)))
  return(nugget)
}
get_cost <- function (param, enugget, esigma2) {
  beta0 <- param[1]
  sigma2 <- exp(param[2])
  x <- rnorm(10^5, beta0, sd = sigma2^0.5)
  prop_bi <- mean(pnorm(x))
  sigma2_bi <- var(pnorm(x))
  sill_bi <- prop * (1-prop)
  nugget_bi <- sill_bi - sigma2_bi
  # return((nugget - enugget) ^ 2 + (esigma2 - sill + nugget) ^ 2)
}
par <- optim(c(-0.2, log(0.1)), get_cost, enugget = 0.18, esigma2 = 0.04)
c(par$par[1], exp(par$par[2]))

get_sill(par$par[1], exp(par$par[2])) - get_nugget(par$par[1], exp(par$par[2]))
get_nugget(par$par[1], exp(par$par[2]))

mesh <- expand.grid(beta0 = seq(-0.7, 0, 0.05), sigma2 = seq(0.5, 1.5, 0.02))
mesh$cost <- purrr::map2_dbl(mesh$beta0, mesh$sigma2,
                             ~ get_cost(c(.x, log(.y)), 0.18, 0.04))

ggplot(mesh, aes(beta0, sigma2)) +
  geom_raster(aes(fill = cost)) +
  scale_fill_viridis_c()


cost2 <- function (param, vg) {
  beta0 <- param[1]
  sigma2 <- exp(param[2])
  phi <- exp(param[3])
  x <- rnorm(10^4, beta0, sd = sigma2^0.5)
  prop_bi <- mean(pnorm(x))
  sigma2_bi <- var(pnorm(x))
  sill_bi <- prop_bi * (1-prop_bi)
  nugget_bi <- sill_bi - sigma2_bi
  theor <- nugget_bi + sigma2_bi * (1 - exp(-vg$dist/phi))
  log(sum(vg$np * (theor - vg$gamma) ^ 2))
  # sum((theor - vg$gamma) ^ 2)
  # theor
}
par <- optim(c(qnorm(mean(data$response)), log(0.03), log(0.03)),
             cost2, vg = vg)
c(par$par[1], exp(par$par[2]), exp(par$par[3]))
# exp(par$par)

cost2(c(-0.5, 0, log(0.02)), vg)

mesh <- expand.grid(beta0 = seq(0, 0.7, 0.01), sigma2 = seq(0.0, 1., 0.01))
mesh$cost <- purrr::map2_dbl(mesh$beta0, mesh$sigma2,
                             ~ cost2(c(.x, log(.y), log(0.02)), vg))
ggplot(mesh, aes(beta0, sigma2)) +
  geom_raster(aes(fill = cost)) +
  # scale_fill_viridis_c(option = "D")
  scale_fill_distiller(palette = "RdYlBu")


plot(vg$gamma, vg0$gamma, asp = 1)
abline(-1.962, 11.801, col =2)

vg_mod <- lm(vg0$gamma ~ vg$gamma)
vg_mod
plot(predict(vg_mod))
lines(vg0$gamma, col = 2)

sd_seq <- seq(0.5, 1.2, 0.01)
plot(sd_seq, sapply(sd_seq, get_variance))

get_prop <- function (c, sigma2) {
  (mean(pnorm(c + rnorm(10^3, sd = sigma2^0.5))) - prop)^2
}
mesh <- expand.grid(c = seq(-1, 0, 0.02), sigma2 = seq(0.5, 1.5, 0.02))
mesh$prop <- purrr::map2_dbl(mesh$c, mesh$sigma2, ~ get_prop(.x, .y))

ggplot(mesh, aes(c, sigma2)) +
  geom_raster(aes(fill = prop)) +
  scale_fill_viridis_c()


cov(
    pnorm(c + rnorm(n, sd = 0.8^2)),
    pnorm(c + rnorm(n, sd = 0.8^2))
    )

cor(
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1),
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1)
    )

table(
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1),
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1)
    )

polycor::polychor(
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1),
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1)
    )

cov(
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1),
    rbinom(n, pnorm(2 + rnorm(n, sd = 0.8^2)), size = 1)
    )

