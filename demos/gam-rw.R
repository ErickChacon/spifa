
# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())
setwd("../src")
getwd()

library(tidyverse)
library(mgcv)
library(raster)
library(RColorBrewer)

Rcpp::sourceCpp("gam.cpp")

# RANDOM WALK 1 ON THE LINE-----------------------------------------------------

n_knots <- 100
D <- diff(diag(n_knots), differences = 1)
tao <- 1
Q <- crossprod(D) / tao ^ 2
rw1 <- rimvnorm_Q_eig(Q)
data.frame(i = 1:n_knots, x = rw1) %>%
  ggplot(aes(i, x)) +
  geom_line() +
  geom_point(col = 2)
mean(rw1)

# RANDOM WALK 2 ON THE LINE ----------------------------------------------------

n_knots <- 100
D <- diff(diag(n_knots), differences = 2)
tao <- 1
Q <- crossprod(D) / tao ^ 2
rw2 <- rimvnorm_Q_eig(Q)
data.frame(i = 1:n_knots, x = rw2) %>%
  ggplot(aes(i, x)) +
  geom_line() +
  geom_point(col = 2)
mean(rw2)

# RANDOM WALK 1 ON A LATTICE ---------------------------------------------------

tao <- 1
n <- 30
r1 <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
ncell1 <- ncell(r1)
adj <- adjacent(r1, cells = 1:ncell1, 4, sorted = TRUE, pairs = TRUE)
adj_tab <- as.numeric(table(adj[, 1]))
Adj <- matrix(0, ncell1, ncell1)
Adj[adj] <- 1
Q <- (diag(adj_tab) - Adj) / tao ^ 2
r1 <- setValues(r1, rimvnorm_Q_eig(Q))
opar <- par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(r1, col = colorRampPalette(brewer.pal(11, "Spectral"))(30), axes = FALSE, box = FALSE)
plot(rasterToPolygons(r1), add = TRUE, border = "black", lwd = 1)
par(opar)

# RANDOM WALK 2 ON A LATTICE ---------------------------------------------------

tao <- 1
n <- 30
r1 <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
ncell1 <- ncell(r1)
index <- 1:ncell1
adj <- adjacent(r1, cells = index, 4, sorted = TRUE, pairs = TRUE)
adj_tab <- as.numeric(table(adj[, 1]))
Adj <- matrix(0, ncell1, ncell1)
Adj[adj] <- 1
diag(Adj) <- - adj_tab
Q <- crossprod(Adj) / tao ^ 2
r1 <- setValues(r1, rimvnorm_Q_eig(Q))
opar <- par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(r1, col = colorRampPalette(brewer.pal(11, "Spectral"))(30), axes = FALSE, box = FALSE)
plot(rasterToPolygons(r1), add = TRUE, border = "black", lwd = 1)
par(opar)



