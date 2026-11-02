#!/usr/bin/env Rscript

# ============================================================================
# Parameter recovery simulation study for stanSPsurv
#
# Runs R replications of each of the 6 model specifications
# (3 structures × 2 baselines) on simulated data with known parameters.
# Reports: bias, RMSE, 90% CI coverage, ESS/sec, divergences, max Rhat.
#
# Usage:
#   Rscript validation/simulation_study.R [--reps 100] [--N 500] [--seed 2026]
# ============================================================================

# -- Setup -------------------------------------------------------------------

local_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

cmdstan_dir <- file.path(getwd(), ".cmdstan", "cmdstan-2.39.0")
if (dir.exists(cmdstan_dir)) cmdstanr::set_cmdstan_path(cmdstan_dir)

library(cmdstanr)

# -- Parse command-line arguments --------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  as.numeric(args[idx + 1L])
}
R     <- get_arg("--reps", 100)
N_obs <- get_arg("--N", 500)
SEED  <- get_arg("--seed", 2026)

message(sprintf("Simulation study: R=%d replications, N=%d observations, seed=%d",
                R, N_obs, SEED))

set.seed(SEED)

# -- True parameter values ---------------------------------------------------

TRUE_BETA  <- c(-1.1, 0.45)
TRUE_GAMMA <- c(0.7, -0.5)
TRUE_SHAPE <- 1.35
TRUE_SIGMA <- 0.35     # frailty SD for exchangeable/spatial

PARAM_NAMES <- c("beta[1]", "beta[2]", "gamma[1]", "gamma[2]", "shape")
TRUE_VALUES <- c(TRUE_BETA, TRUE_GAMMA, TRUE_SHAPE)

# -- DGP ---------------------------------------------------------------------

simulate_data <- function(structure, baseline, N = N_obs, J = 8L) {
  x <- rnorm(N)
  z <- rnorm(N)
  X <- cbind(1, x)
  Z <- cbind(1, z)

  group_id <- rep(seq_len(J), length.out = N)

  # Frailties (survival stage only, for simplicity)
  u <- rep(0, J)
  if (structure == "exchangeable") {
    u <- rnorm(J, 0, TRUE_SIGMA)
    u <- u - mean(u)
  } else if (structure == "spatial") {
    # Spatially correlated frailties on a chain graph
    u[1] <- rnorm(1, 0, TRUE_SIGMA)
    for (j in 2:J) u[j] <- 0.6 * u[j - 1] + rnorm(1, 0, TRUE_SIGMA * 0.5)
    u <- u - mean(u)
  }

  pi <- plogis(drop(Z %*% TRUE_GAMMA))
  susceptible <- rbinom(N, 1, pi)

  eta <- drop(X %*% TRUE_BETA) + u[group_id]
  q <- runif(N)
  if (baseline == "weibull") {
    latent <- (-log1p(-q) / exp(eta))^(1 / TRUE_SHAPE)
  } else {
    latent <- (q / ((1 - q) * exp(eta)))^(1 / TRUE_SHAPE)
  }
  latent[susceptible == 0] <- Inf
  censor <- rexp(N, rate = 0.10)
  time <- pmin(latent, censor)
  event <- as.integer(latent <= censor)

  # Adjacency: chain graph

  A <- matrix(0, J, J)
  for (j in seq_len(J - 1L)) { A[j, j + 1L] <- 1; A[j + 1L, j] <- 1 }

  edges <- which(upper.tri(A) & A == 1, arr.ind = TRUE)
  graph <- stanSPsurv:::adjacency_graph(if (structure == "spatial") A else NULL)

  stan_data <- list(
    N = N, M = N, K = ncol(X), Q = ncol(Z),
    t0 = rep(0, N), time = time, event = as.integer(event),
    X = X, Z_spell = Z,
    spell_start = seq_len(N), spell_end = seq_len(N),
    spell_event = as.integer(event),
    J = if (structure == "pooled") 1L else J,
    spell_group_id = if (structure == "pooled") rep(1L, N) else group_id,
    N_edges = if (structure == "spatial") nrow(edges) else 0L,
    node1 = if (structure == "spatial") as.integer(edges[, 1]) else integer(),
    node2 = if (structure == "spatial") as.integer(edges[, 2]) else integer(),
    N_components = graph$N_components,
    component_id = graph$component_id,
    component_size = graph$component_size,
    icar_scale = graph$scale,
    baseline = match(baseline, c("weibull", "loglogistic"))
  )

  list(stan_data = stan_data, event_rate = mean(event))
}

# -- Compile Stan models once ------------------------------------------------

message("Compiling Stan models...")
models <- lapply(c("pooled", "exchangeable", "spatial"), function(s) {
  cmdstan_model(file.path("inst", "stan", paste0(s, "_spsurv.stan")),
                force_recompile = FALSE, quiet = TRUE)
})
names(models) <- c("pooled", "exchangeable", "spatial")
message("Compilation complete.")

# -- Experimental conditions -------------------------------------------------

conditions <- expand.grid(
  structure = c("pooled", "exchangeable", "spatial"),
  baseline  = c("weibull", "loglogistic"),
  stringsAsFactors = FALSE
)

# -- Run simulations ---------------------------------------------------------

out_dir <- file.path("validation", "simulation_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- vector("list", nrow(conditions) * R)
idx <- 0L

for (cond in seq_len(nrow(conditions))) {
  structure <- conditions$structure[cond]
  baseline  <- conditions$baseline[cond]
  message(sprintf("\n=== %s / %s ===", structure, baseline))

  for (r in seq_len(R)) {
    idx <- idx + 1L
    if (r %% 10 == 0 || r == 1) message(sprintf("  replication %d/%d", r, R))

    sim <- simulate_data(structure, baseline)

    t0 <- proc.time()
    fit <- tryCatch(
      models[[structure]]$sample(
        data = sim$stan_data,
        seed = SEED * 1000 + idx,
        chains = 4, parallel_chains = 4,
        iter_warmup = 500, iter_sampling = 500,
        refresh = 0,
        adapt_delta = 0.90, max_treedepth = 11,
        show_messages = FALSE, show_exceptions = FALSE
      ),
      error = function(e) NULL
    )
    elapsed <- (proc.time() - t0)[["elapsed"]]

    if (is.null(fit)) {
      all_results[[idx]] <- data.frame(
        structure = structure, baseline = baseline, rep = r,
        param = PARAM_NAMES, truth = TRUE_VALUES,
        mean = NA, q05 = NA, q95 = NA,
        ess_bulk = NA, rhat = NA,
        divergences = NA, max_treedepth_hits = NA,
        elapsed_sec = elapsed, event_rate = sim$event_rate,
        stringsAsFactors = FALSE
      )
      next
    }

    draws <- fit$draws(variables = c("beta", "gamma", "shape"),
                       format = "matrix")
    diag <- fit$diagnostic_summary(quiet = TRUE)
    divs <- sum(diag$num_divergent)
    mtd  <- sum(diag$num_max_treedepth)

    sel <- match(PARAM_NAMES, colnames(draws))
    d <- draws[, sel, drop = FALSE]
    est_mean <- colMeans(d)
    est_q05  <- apply(d, 2, quantile, 0.05)
    est_q95  <- apply(d, 2, quantile, 0.95)
    est_ess  <- apply(d, 2, posterior::ess_bulk)
    est_rhat <- apply(d, 2, posterior::rhat)

    all_results[[idx]] <- data.frame(
      structure = structure, baseline = baseline, rep = r,
      param = PARAM_NAMES, truth = TRUE_VALUES,
      mean = est_mean, q05 = est_q05, q95 = est_q95,
      ess_bulk = est_ess, rhat = est_rhat,
      divergences = divs, max_treedepth_hits = mtd,
      elapsed_sec = elapsed, event_rate = sim$event_rate,
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, all_results)
write.csv(results, file.path(out_dir, "raw_results.csv"), row.names = FALSE)

# -- Summary table -----------------------------------------------------------

message("\nComputing summary statistics...")

summary_table <- do.call(rbind, lapply(
  split(results, interaction(results$structure, results$baseline, results$param)),
  function(df) {
    df <- df[!is.na(df$mean), ]
    if (nrow(df) == 0) return(NULL)
    bias     <- mean(df$mean - df$truth)
    rmse     <- sqrt(mean((df$mean - df$truth)^2))
    coverage <- mean(df$q05 <= df$truth & df$truth <= df$q95)
    med_ess  <- median(df$ess_bulk)
    med_rhat <- median(df$rhat)
    med_time <- median(df$elapsed_sec)
    ess_sec  <- med_ess / med_time
    tot_divs <- sum(df$divergences)

    data.frame(
      structure = df$structure[1], baseline = df$baseline[1],
      param = df$param[1], truth = df$truth[1],
      bias = bias, rmse = rmse, coverage_90 = coverage,
      median_ess = med_ess, median_rhat = med_rhat,
      ess_per_sec = ess_sec, total_divergences = tot_divs,
      n_converged = nrow(df), median_seconds = med_time,
      stringsAsFactors = FALSE
    )
  }
))
rownames(summary_table) <- NULL
write.csv(summary_table, file.path(out_dir, "summary.csv"), row.names = FALSE)

message("\n=== SUMMARY ===")
print(summary_table, digits = 3, row.names = FALSE)
message(sprintf("\nResults saved to %s/", out_dir))
