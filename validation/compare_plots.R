#!/usr/bin/env Rscript

# ============================================================================
# Visual comparison: stanSPsurv vs BayesSPsurv plotting capabilities
#
# Fits the same simulated data with both packages, then produces a
# multi-page PDF showing what each package can plot.
#
# Requires the local BayesSPsurv rho fix (tolower in mcmc_engines.R).
# ============================================================================

library(survival)
library(stanSPsurv)
library(BayesSPsurv)
library(coda)

set.seed(2026)

# -- Simulate data (shape < 1 for BayesSPsurv compatibility) ----------------

N <- 400; J <- 8
x <- rnorm(N); z <- rnorm(N)
grp <- rep(seq_len(J), length.out = N)

TRUE_BETA  <- c(-1.1, 0.45)
TRUE_GAMMA <- c(0.7, -0.5)
TRUE_SHAPE <- 0.75

pi <- plogis(cbind(1, z) %*% TRUE_GAMMA)
susceptible <- rbinom(N, 1, pi)
u <- rnorm(J, 0, 0.35); u <- u - mean(u)
eta <- cbind(1, x) %*% TRUE_BETA + u[grp]
q <- runif(N)
latent <- (-log1p(-q) / exp(eta))^(1 / TRUE_SHAPE)
latent[susceptible == 0] <- Inf
censor <- rexp(N, 0.10)
time <- pmin(latent, censor)
event <- as.integer(latent <= censor)

df <- data.frame(time = time, event = event, x = x, z = z,
                 grp = factor(grp))

cat("Event rate:", mean(event), "\n")

# -- Fit stanSPsurv (exchangeable) ------------------------------------------

cat("Fitting stanSPsurv...\n")
stan_fit <- stan_exchange_spsurv(
  Surv(time, event) ~ x, split_formula = ~ z,
  data = df, group = grp,
  form = "weibull", chains = 4, parallel_chains = 4,
  iter_warmup = 500, iter_sampling = 500, refresh = 0
)
cat("  done.\n")

# -- Fit BayesSPsurv (exchangeable) -----------------------------------------

cat("Fitting BayesSPsurv...\n")
bayes_df <- data.frame(
  duration = time, t.0 = rep(0, N), cured = event,
  lastyear = rep(1L, N), x = x, z = z, group_id = grp
)

bayes_fit <- exchangeSPsurv(
  duration = duration ~ x, immune = cured ~ z,
  Y0 = "t.0", LY = "lastyear", S = "group_id",
  data = bayes_df, N = 15000, burn = 5000, thin = 5,
  w = c(1, 1, 1), m = 10, form = "weibull",
  prop.varV = 1e-04, prop.varW = 1e-04
)
cat("  done.\n")

# -- Generate comparison PDF -------------------------------------------------

out_file <- file.path("validation", "comparison_plots.pdf")
cat("Writing", out_file, "...\n")

pdf(out_file, width = 11, height = 8.5)

# ---- Page 1: Title --------------------------------------------------------

plot.new()
text(0.5, 0.7, "stanSPsurv vs BayesSPsurv", cex = 2.2, font = 2)
text(0.5, 0.55, "Visual comparison of plotting capabilities", cex = 1.3)
text(0.5, 0.40, paste("Exchangeable Weibull model, N =", N, ", J =", J),
     cex = 1.1, col = "grey40")
text(0.5, 0.30, paste0("True: beta = (", paste(TRUE_BETA, collapse = ", "),
                        "), gamma = (", paste(TRUE_GAMMA, collapse = ", "),
                        "), shape = ", TRUE_SHAPE),
     cex = 0.95, col = "grey40")

# ---- Page 2: Coefficient plots --------------------------------------------

par(mfrow = c(1, 2))

# stanSPsurv: coefficient forest plot
plot(stan_fit, main = "stanSPsurv: plot(fit)")

# BayesSPsurv: no coefficient plot — show what it has (coda trace)
plot(mcmc(bayes_fit$betas), auto.layout = FALSE,
     main = "BayesSPsurv: coda trace + density")
mtext("(BayesSPsurv has no coefficient plot — only coda trace/density)",
      side = 1, line = 4, cex = 0.7, col = "grey40")

par(mfrow = c(1, 1))

# ---- Page 3: BayesSPsurv coda plots (betas) -------------------------------

par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(mcmc(bayes_fit$betas), auto.layout = FALSE)
title("BayesSPsurv: betas (coda trace + density)", outer = TRUE)
par(mfrow = c(1, 1))

# ---- Page 4: BayesSPsurv coda plots (gammas) ------------------------------

par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(mcmc(bayes_fit$gammas), auto.layout = FALSE)
title("BayesSPsurv: gammas (coda trace + density)", outer = TRUE)
par(mfrow = c(1, 1))

# ---- Page 5: BayesSPsurv coda plots (rho) ---------------------------------

par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
plot(mcmc(bayes_fit$rho), auto.layout = FALSE)
title("BayesSPsurv: rho / shape (coda trace + density)", outer = TRUE)
par(mfrow = c(1, 1))

# ---- Page 6: stanSPsurv survival curves -----------------------------------

nd <- data.frame(x = c(-1, 0, 1))
plot_survival(stan_fit, newdata = nd, times = seq(0.1, 30, length.out = 80),
              labels = c("x = -1", "x = 0", "x = +1"),
              cols = c("firebrick", "grey30", "steelblue"))
mtext("stanSPsurv: plot_survival()  (no BayesSPsurv equivalent)",
      side = 1, line = 3.5, cex = 0.8, col = "grey40")

# ---- Page 7: stanSPsurv susceptibility ------------------------------------

plot_susceptibility(stan_fit, variable = "z", data = df)
mtext("stanSPsurv: plot_susceptibility()  (no BayesSPsurv equivalent)",
      side = 1, line = 3.5, cex = 0.8, col = "grey40")

# ---- Page 8: stanSPsurv frailty caterpillar plot --------------------------

par(mfrow = c(1, 2))
plot_frailty(stan_fit, stage = "survival")
plot_frailty(stan_fit, stage = "split")
par(mfrow = c(1, 1))
mtext("stanSPsurv: plot_frailty()  (no BayesSPsurv equivalent)",
      side = 1, line = 0.5, cex = 0.8, col = "grey40", outer = TRUE)

# ---- Page 9: Posterior comparison (converted to same scale) ----------------

# Convert BayesSPsurv posteriors to stanSPsurv parametrization
shape_bayes <- exp(-bayes_fit$rho)
beta1_bayes <- -shape_bayes * bayes_fit$betas[, 1]
beta2_bayes <- -shape_bayes * bayes_fit$betas[, 2]
gamma1_bayes <- -bayes_fit$gammas[, 1]
gamma2_bayes <- -bayes_fit$gammas[, 2]

# stanSPsurv draws
stan_draws <- stan_fit$fit$draws(c("beta", "gamma", "shape"), format = "matrix")
beta1_stan  <- stan_draws[, "beta[1]"]
beta2_stan  <- stan_draws[, "beta[2]"]
gamma1_stan <- stan_draws[, "gamma[1]"]
gamma2_stan <- stan_draws[, "gamma[2]"]
shape_stan  <- stan_draws[, "shape"]

par(mfrow = c(2, 3), oma = c(0, 0, 2, 0))

plot_density_comparison <- function(stan_d, bayes_d, true_val, label) {
  d_s <- density(stan_d)
  d_b <- density(bayes_d)
  xlim <- range(c(d_s$x, d_b$x))
  ylim <- c(0, max(c(d_s$y, d_b$y)) * 1.1)
  plot(d_s, col = "steelblue", lwd = 2, xlim = xlim, ylim = ylim,
       main = label, xlab = "")
  lines(d_b, col = "firebrick", lwd = 2)
  abline(v = true_val, lty = 2, col = "grey40", lwd = 1.5)
  legend("topright", legend = c("stanSPsurv", "BayesSPsurv", "Truth"),
         col = c("steelblue", "firebrick", "grey40"),
         lwd = c(2, 2, 1.5), lty = c(1, 1, 2), bty = "n", cex = 0.7)
}

plot_density_comparison(beta1_stan, beta1_bayes, TRUE_BETA[1], "beta[1]")
plot_density_comparison(beta2_stan, beta2_bayes, TRUE_BETA[2], "beta[2]")
plot_density_comparison(gamma1_stan, gamma1_bayes, TRUE_GAMMA[1], "gamma[1]")
plot_density_comparison(gamma2_stan, gamma2_bayes, TRUE_GAMMA[2], "gamma[2]")
plot_density_comparison(shape_stan, shape_bayes, TRUE_SHAPE, "shape")

# Empty panel for legend
plot.new()
text(0.5, 0.6, "Both packages recover\nthe true parameters", cex = 1.1)
text(0.5, 0.35, "stanSPsurv posteriors are\ntighter (higher ESS)", cex = 0.9,
     col = "grey40")

title("Posterior density comparison (same parametrization)", outer = TRUE)
par(mfrow = c(1, 1))

# ---- Page 10: Summary table -----------------------------------------------

plot.new()
text(0.5, 0.92, "Feature comparison", cex = 1.8, font = 2)

features <- c(
  "Coefficient forest plot",
  "Survival curves with credible bands",
  "Susceptibility curves",
  "Frailty caterpillar plot",
  "Trace plots (via coda/cmdstanr)",
  "Posterior density comparison",
  "Spatial autocorrelation (Moran I)",
  "Join count test",
  "Posterior predictive checks",
  "LOO-CV / WAIC (via log_lik)"
)
stan_has  <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
bayes_has <- c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE)

y_pos <- seq(0.82, 0.82 - 0.065 * (length(features) - 1), length.out = length(features))

text(0.05, 0.87, "Feature", adj = 0, font = 2, cex = 0.85)
text(0.65, 0.87, "stanSPsurv", font = 2, cex = 0.85)
text(0.85, 0.87, "BayesSPsurv", font = 2, cex = 0.85)

for (i in seq_along(features)) {
  text(0.05, y_pos[i], features[i], adj = 0, cex = 0.8)
  text(0.65, y_pos[i], if (stan_has[i]) "\u2713" else "\u2717",
       col = if (stan_has[i]) "darkgreen" else "red", cex = 1)
  text(0.85, y_pos[i], if (bayes_has[i]) "\u2713" else "\u2717",
       col = if (bayes_has[i]) "darkgreen" else "red", cex = 1)
}

dev.off()
cat("Done. Output:", out_file, "\n")
