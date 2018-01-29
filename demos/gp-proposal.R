
library(datasim)
library(spmirt)

# SIMULATE GAUSSIAN PROCESS ----------------------------------------------------

n <- 300
f <- list(
  mean ~ I(0) +
    gp(list(s1, s2), cor.model = "exp_cor", cor.params = list(phi = 0.02),
       sigma2 = 0.5),
  sd ~ I(0)
  )
data <- sim_model(formula = f,
                  link_inv = list(identity, identity),
                  generator = rnorm,
                  n = n)

# FUNCTION TO BE EVALUATED -----------------------------------------------------

dist <- as.matrix(dist(dplyr::select(data, s1, s2)))
Sigma <- 2 * exp(-dist/0.02)
fun <- function(x) {
  # sigma2 <- exp(x[1])
  # phi <- exp(x[2])
  sigma2 <- psych::logistic(x[1])
  phi <- exp(x[2])
  Sigma <- sigma2 * exp(-dist / phi)
  spmirt::dmvnorm(data$gp.list.mean, mean=rep(0, n), sigma=Sigma)
}

# QUADRATIC APPROXIMATION ------------------------------------------------------

npts <- 30
argRanges <- list(c(-0.5, 0.5), c(-5, -3))
(approx <- spatsurv::QuadApprox(fun, npts, argRanges, plot = TRUE))
(varmat <- solve(-approx$curvature))
#            [,1]       [,2]
# [1,] 0.03276579 0.01591923
# [2,] 0.01591923 0.04475116
D <- diag(1 / diag(varmat) ^ 0.5)
D %*% varmat %*% D
c(psych::logistic(approx$max[1]), exp(approx$max[2]))

# psych::logit(0.5) # [1] 0
# log(0.02) # [1] -3.912023
# # lapack dgesv

vg <- gstat::variogram(gp.list.mean ~ 1, ~ s1 + s2, data, cutoff = 0.4, width = 0.005)
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  # geom_smooth() +
  # geom_smooth(aes(weight = np, fill = 2)) +
  geom_smooth(method = mgcv::gam, formula = y ~ s(x), se = FALSE) +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.4))

