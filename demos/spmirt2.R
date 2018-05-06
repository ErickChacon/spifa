# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())
library(datasim)
library(tidyverse)

# SIMULATE MULTIVARIATE SPATIAL DATA -------------------------------------------

# set.seed(4)
Corr <- matrix(c(1, -0.3, 0, -0.3, 1, 0.3, 0, 0.3, 1), nrow = 3)
sigmas <- rep(0.4^0.5, 3)
D <- diag(sigmas)
Cov <- D %*% Corr %*% D

# beta <- c(-0.5, 0, 0.5)
beta <- c(0, 0, 0)
variance <- 0.6 * matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)
cor.model <- "exp_cor"
cor.params <- list(list(phi = 0.04), list(phi = 0.04), list(phi = 0.1))

f <- list(
  mean ~ mfe(x1, beta = get("beta")) +
    mre(factor(id), sigma = get("Cov")) +
    mgp(list(s1), variance = get("variance"), cor.model = get("cor.model"),
        cor.params = get("cor.params")),
  sd ~ I(0)
  )

n <- 300
m <- 3
(data_geo <- sim_model(formula = f, n = n, responses = m))
# knitr::kable(head(data_model, 10))

X <- matrix(rnorm(20), 10, 2)

# VISUALIZE MULTIVARIATE SPATIAL DATA ------------------------------------------

ggplot(data_geo, aes(x1, response)) +
  geom_smooth(aes(col = factor(response_label))) +
  geom_point(aes(col = factor(response_label)))

ggplot(data_geo, aes(s1, mgp.list.mean)) +
  geom_line(aes(col = factor(response_label)))

data_geo %>%
  dplyr::select(id, mre.factor.mean, response_label) %>%
  spread(response_label, mre.factor.mean) %>%
  dplyr::select(-id) %>%
  GGally::ggpairs(aes(fill = "any"))

data_geo_wide <- data_geo %>%
  dplyr::rename(ability = response, id_person = id) %>%
  gather(var, value, mre.factor.mean:ability) %>%
  mutate(var = paste0(var, response_label)) %>%
  select(-response_label) %>%
  spread(var, value)


# SIMULATE ITEM FACTOR DATA ----------------------------------------------------

q <- 10
init_data <- purrr::map(1:q, ~ data_geo_wide) %>%
  purrr::reduce(rbind)

# n <- 300
difficulty <- matrix((1:q - 5)/10 * 2, nrow = 1)
discrimination1 <- seq(0.4, 1.5, length.out = q)
discrimination2 <- runif(q, 0, 2)
discrimination3 <- runif(q, 0, 2)
discrimination1[1] <- 1
discrimination1[c(3, 5, 8)] <- 0
discrimination2[1:2] <- c(0, 1)
discrimination2[c(4, 5, 10)] <- 0
# discrimination3[1:3] <- c(0, 0, 1)
# discrimination1 <- discrimination1 * 0.3
# discrimination2 <- discrimination2 * 0.3
cbind(discrimination1, discrimination2, discrimination3)

f <- list(
  prob ~ mfa(ones, beta = get("difficulty")) +
    mfe(ability1, beta = get("discrimination1")) +
    mfe(ability2, beta = get("discrimination2")),
  # + mfe(ability3, beta = get("discrimination3")),
  size ~ I(1)
  )

data_long <- sim_model(formula = f,
                        link_inv = list(pnorm, identity),
                        generator = rbinom,
                        responses = q,
                        n = n,
                        init_data = init_data
                        )

data_long <- dplyr::rename(data_long, subject = id,
                           item = response_label, y = response)

# VISUALIZE ITEM FACTOR DATA ---------------------------------------------------

explor <- data_long %>%
  group_by(subject) %>%
  summarize(endorse = mean(y),
            ability1 = unique(ability1),
            ability2 = unique(ability2),
            # ability3 = unique(ability3),
            x1 = unique(x1))
ggplot(explor, aes(ability1, endorse)) + geom_point(alpha = 0.5)
ggplot(explor, aes(ability2, endorse)) + geom_point(alpha = 0.5)
# ggplot(explor, aes(ability3, endorse)) + geom_point(alpha = 0.5)
# ggplot(explor, aes(x1, endorse)) + geom_point(alpha = 0.5)

# PREPARE DATA -----------------------------------------------------------------

response <- data_long$y
dist <- as.matrix(dist(dplyr::select(data_geo_wide, s1)))
# dist <- as.matrix(dist(dplyr::select(data_geo_wide, s1)[order(data_geo_wide$s1),]))
# dist <- dist[order(data_geo_wide$s1),]
n
q
m <- 2
iter <- 5 * 10 ^ 4
# iter <- 5 * 10 ^ 2
cor.params <- c(0.04, 0.04)
sig.params <- c(0.6 ^ 0.5, 0.6 ^ 0.5)
fix.sigma <- 0.4^0.5
# sigma_prop <- matrix(c(0.138, -0.023, -0.023, 0.1), 2) * 2.38 ^ 2 / 2
sigma_prop <- 0.001 * diag(5)
disc_mat <- cbind(discrimination1, discrimination2)
L_a <- lower.tri(disc_mat, diag = TRUE) * 1
T_gp <- diag(m)

# RUN --------------------------------------------------------------------------

Rcpp::sourceCpp("../src/mirt-gibss-sp.cpp")
source("../R/ggplot-mcmc.R")
Rcpp::sourceCpp("../src/ifa-main.cpp")

# set.seed(5)
system.time(
  samples <- ifa_gibbs_sp(response, dist, n, q, m, cor.params, sig.params,
                          Corr[1:2, 1:2], fix.sigma, sigma_prop, L_a, T_gp, 0.234,
                          iter)
)

# system.time(
#   samples <- spmirt(response = response, nobs = n, nitems = q, nfactors = 2,
#                     L_rest = L_a, niter = iter)
#   )

samples_tib <- as_tibble.spmirt.list(samples, iter/2)
summary(samples_tib)
samples_long <- gather(samples_tib)

as_tibble.spmirt.list(samples, 0, 10, "c") %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples, 0, 10, "a") %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples, iter/2, 10, "a") %>%
  gg_density(alpha = 0.5, ridges = TRUE, aes(fill = Parameters), scale = 4)

as_tibble.spmirt.list(samples, iter/2, 10, "theta") %>%
  dplyr::select(1:100) %>%
  gg_density(alpha = 0.5, ridges = TRUE, aes(fill = Parameters), scale = 4)

as_tibble.spmirt.list(samples, 0, 10, "theta") %>%
  select(1:10) %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples, 0, 10, "mgp_sd") %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples, 0, 10, "mgp_phi") %>%
  gg_trace(alpha = 0.6)

as_tibble.spmirt.list(samples, 0, 10, "a") %>%
  gg_density2d(`Discrimination 1`, `Discrimination 2`, each = 10,
               keys = c("Item ", "Discrimination "),
               highlight = c(discrimination1, discrimination2))

as_tibble.spmirt.list(samples, 0, 10, "a") %>%
  gg_scatter(`Discrimination 1`, `Discrimination 2`, each = 10,
               keys = c("Item ", "Discrimination "),
               highlight = c(discrimination1, discrimination2))

as_tibble.spmirt.list(samples, iter/ 2, select = "a") %>%
  summary() %>%
  mutate(param = c(discrimination1, discrimination2)) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)

as_tibble.spmirt.list(samples, iter/2, select = "c") %>%
  summary() %>%
  mutate(param = as.numeric(difficulty)) %>%
  gg_errorbarh() +
  geom_point(aes(param, Parameters), col = 3)

as_tibble.spmirt.list(samples, iter/2, select = "theta") %>%
  dplyr::select(1:300) %>%
  summary() %>%
  mutate(param = data_geo$response[1:300]) %>%
  gg_errorbarh(sorted = TRUE) +
  geom_point(aes(x = param), col = 3)

as_tibble.spmirt.list(samples, iter/2, select = "theta") %>%
  dplyr::select(301:600) %>%
  summary() %>%
  mutate(param = data_geo$response[301:600]) %>%
  gg_errorbarh(sorted = TRUE) +
  geom_point(aes(x = param), col = 3)

ability1_pred <- as_tibble.spmirt.list(samples, iter/2, select = "theta") %>%
  dplyr::select(1:300) %>%
  summary() %>%
  mutate(param = data_geo$response[1:300],
         s1 = data_geo$s1[1:300],
         s2 = s1,
         estim = `50%`)
ability1_pred %>%
    ggplot(aes(s1, `50%`)) +
    geom_line() +
    geom_line(aes(s1, param, col = "real"))

vg <- gstat::variogram(estim ~ 1, ~ s1 + s2, ability1_pred, cutoff = 1, width = 0.01)
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.7))

ability2_pred <- as_tibble.spmirt.list(samples, iter/2, select = "theta") %>%
  dplyr::select(301:600) %>%
  summary() %>%
  mutate(param = data_geo$response[301:600],
         s1 = data_geo$s1[301:600],
         s2 = s1,
         estim = `50%`)
ability2_pred %>%
  ggplot(aes(s1, `50%`)) +
  geom_line() +
  geom_line(aes(s1, param, col = "real"))

vg <- gstat::variogram(estim ~ 1, ~ s1 + s2, ability2_pred, cutoff = 1, width = 0.005)
ggplot(vg, aes(dist, gamma)) +
  geom_point(aes(size = np)) +
  geom_smooth() +
  expand_limits(y = 0, x = 0) +
  scale_x_continuous(limits = c(0, 0.7))

# # # PREPARE DATA FOR MODELLING ---------------------------------------------------
# #
# # Y <- data_model %>% dplyr::select(id, response, response_label) %>%
# #   spread(response_label, response) %>%
# #   arrange(id) %>%
# #   dplyr::select(-id) %>%
# #   as.matrix()
# #
# # X <- data_model %>% dplyr::select(id, matches("^x[[:digit:]]+$")) %>%
# #   unique() %>%
# #   arrange(id) %>%
# #   dplyr::select(-id) %>%
# #   as.matrix()
# #
# # Beta <- matrix(beta, nrow = 1)
# # Sigma_proposal <- diag(1, 3)
# #
# # # RUN MODEL --------------------------------------------------------------------
# #
# # getwd()
# # Rcpp::sourceCpp("../src/multi-lm.cpp")
# # source("../R/ggplot-mcmc.R")
# #
# # iter <- 10^6
# # system.time(
# #   samples <- multi_lm(Y, X, iter, 0.01 * Sigma_proposal, 0.001 * Sigma_proposal)
# # )
# # samples %>% map(~ tail(.))
# #
# # # apply(samples$beta, 2, mean)
# # # cor(samples$beta)
# #
# # # Visualize traces
# # as_tibble(samples, 0, 100, select = "beta") %>%
# #   gg_trace(wrap = TRUE, alpha = 0.6)
# #
# # as_tibble(samples, 0, 100, select = "beta") %>% gg_trace(alpha = 0.6)
# # as_tibble(samples, 0, 100, select = "corr_chol") %>% gg_trace(alpha = 0.6)
# # as_tibble(samples, 0, 100, select = "corr") %>% gg_trace(alpha = 0.6)
# # as_tibble(samples, 0, 100, select = "sigmas") %>% gg_trace(alpha = 0.6)
# #
# # bla <- as_tibble(samples, iter/2, select = "sigmas")
# # cov(log(bla))
# # nrow(unique(bla)) / nrow(bla)
# #
# # bla <- as_tibble(samples, iter/2, select = "corr_chol")
# # cov(bla)
# # nrow(unique(bla)) / nrow(bla)
# #
# # # Visualize densities
# #
# # as_tibble(samples, iter / 2, select = "corr_chol") %>%
# #   gg_density(aes(fill = Parameters), scale = 2, alpha = 0.5, ridges = TRUE)
# #
# # as_tibble(samples, iter / 2, select = "corr") %>%
# #   gg_density(aes(fill = Parameters), scale = 1, alpha = 0.5, ridges = TRUE)
# #
# # # Visualize credible intervals
# # as_tibble(samples, iter / 2, select = "beta") %>%
# #   summary() %>%
# #   mutate(param = beta) %>%
# #   gg_errorbarh() +
# #   geom_point(aes(param, Parameters), col = 3)
# #
# # Corr_chol <- t(chol(Corr))
# # corr_chol <- Corr_chol[lower.tri(Corr_chol, diag = TRUE)]
# # corr <- Corr[lower.tri(Corr)]
# #
# # as_tibble(samples, iter / 2, select = "corr_chol") %>%
# #   summary() %>%
# #   mutate(param = corr_chol) %>%
# #   gg_errorbarh() +
# #   geom_point(aes(param, Parameters), col = 3)
# #
# # as_tibble(samples, iter / 2, select = "corr") %>%
# #   summary() %>%
# #   mutate(param = corr) %>%
# #   gg_errorbarh() +
# #   geom_point(aes(param, Parameters), col = 3)
# #
# #
# # as_tibble(samples, iter / 2 ,select = "sigmas") %>%
# #   summary() %>%
# #   mutate(param = sigmas) %>%
# #   gg_errorbarh() +
# #   geom_point(aes(param, Parameters), col = 3)
# #
# #
# # # Visualize credible intervals for all Parameters
# # as_tibble(samples, iter / 2) %>%
# #   summary() %>%
# #   mutate(param = c(beta, corr_chol, corr, sigmas)) %>%
# #   gg_errorbar() +
# #   geom_point(aes(Parameters, param), col = 3)
# #
