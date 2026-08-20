# Simulate the `ipixuna` dataset shipped with the package: simulated household
# food-insecurity survey data for the urba area of Ipixuna (Para state, Brazil). An sf
# object with an `items` response matrix and a spatial locations. See R/spifa-package.R
# for the shipped documentation.
#
#   theta_ik = effect_k wealth_i + w_k(coords_i) + xi_ik,  xi_i ~ N(0, Cov)
#   w_k(.) ~ GP(0, A_kk^2 * exp(-dist / phi_k))
#   z_ij = easiness_j + a_j'theta_i + eps_ij,  eps_ij ~ N(0, 1)
#   items_ij = 1{z_ij > 0}

## Setup and dimensions

library(sf)
library(dplyr)

set.seed(202608)

n <- 100  # households
q <- 10   # items
nfactors <- 3    # latent factors
ngp <- 3    # independent GPs (one dedicated GP per factor)

## Ifa parameters

easiness <- seq(-0.68, 0.85, length.out = q)

discrimination <- matrix(runif(q * nfactors, 0, 1.5), q)
discrimination[discrimination < 0.5] <- 0
discrimination[, 3] <- - discrimination[, 3]

## Latent abilities

# Predictor effect
wealth <- rnorm(n)
effect <- c(-0.7, 0.7, 0.3)
wealth_effect <- outer(wealth, effect)

# Non-spatial term
resid_params <- list(sd = sqrt(c(0.2, 0.3, 0.2)), corr = c(-0.7, 0.1, -0.1))

nsp_corr <- diag(nfactors)
nsp_corr[lower.tri(nsp_corr)] <- resid_params$corr
nsp_corr[upper.tri(nsp_corr)] <- resid_params$corr
nsp_cov <- diag(resid_params$sd) %*% nsp_corr %*% diag(resid_params$sd)
nsp_chol <- chol(nsp_cov)
nonspatial <- matrix(rnorm(n * nfactors), n) %*% nsp_chol

# Spatial term
mgp_params <- list(sd = sqrt(c(0.4, 0.3, 0.5)), phi = c(100, 150, 250))

bbox <- c(xmin = -71.70038, ymin = -7.06058, xmax = -71.68109, ymax = -7.03724)
# households are concentrated towards the town centre rather than uniform
# across the bounding box
coords <- data.frame(
    lon = bbox["xmin"] + rbeta(n, 4, 4) * (bbox["xmax"] - bbox["xmin"]),
    lat = bbox["ymin"] + rbeta(n, 4, 4) * (bbox["ymax"] - bbox["ymin"])) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_geometry()
distances <- st_distance(coords) |> units::drop_units()
# mean(as.numeric(distances) < 300)
# mean(as.numeric(distances) < 450)
# mean(as.numeric(distances) < 750)

spatial <- sapply(seq_len(ngp), function (k) {
  gp_cov <- mgp_params$sd[k]^2 * exp(-distances / mgp_params$phi[k])
  as.numeric(chol(gp_cov) %*% rnorm(n))
})

# Abilities

abilities <- wealth_effect + nonspatial + spatial

## Items

z <- outer(rep(1, n), easiness) + abilities %*% t(discrimination) + matrix(rnorm(n * q), n)
items <- (z > 0) * 1

## Ipixuna dataset
ipixuna <- st_sf(id = seq_len(n), wealth = wealth, items = I(items), geometry = coords)
attr(ipixuna, "parameters") <- list(
    easiness = easiness, discrimination = discrimination, abilities = abilities,
    effect = effect, resid_params = resid_params, mgp_params = mgp_params
)

save(ipixuna, file = "data/ipixuna.RData")
