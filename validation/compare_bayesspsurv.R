#!/usr/bin/env Rscript

# ============================================================================
# Formal comparison: stanSPsurv vs BayesSPsurv
#
# Fits the same simulated datasets with both packages and compares
# bias, RMSE, 90% CI coverage, effective sample size, and wall time.
#
# Parametrization mapping (derived from the C++ likelihood in BayesSPsurv):
#   beta_stan  = -shape * beta_bayes    (BayesSPsurv uses eXB = exp(-X*beta))
#   gamma_stan = -gamma_bayes           (stanSPsurv pi = P(susceptible);
#                                        BayesSPsurv delta = P(immune))
#   shape_stan = exp(-rho_bayes)
#
# IMPORTANT NOTES:
#   - BayesSPsurv constrains rho >= 0.01 with Gamma(1,1) prior, so
#     shape = exp(-rho) <= 0.99.  Only decreasing hazards (shape < 1)
#     are within its support.  We use shape = 0.75 for a fair comparison.
#   - BayesSPsurv v0.1.5 has a case-mismatch bug: the MCMC engine checks
#     `form %in% "Weibull"` but receives `"weibull"`, so rho is never
#     sampled.  This comparison requires the local fix (tolower).
#
# Usage:
#   Rscript validation/compare_bayesspsurv.R [--reps 50] [--N 400] [--seed 2026]
# ============================================================================

# -- Setup -------------------------------------------------------------------

local_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

cmdstan_dir <- file.path(getwd(), ".cmdstan", "cmdstan-2.39.0")
if (dir.exists(cmdstan_dir)) cmdstanr::set_cmdstan_path(cmdstan_dir)

library(cmdstanr)
library(BayesSPsurv)
library(coda)

# -- Parse command-line arguments --------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  as.numeric(args[idx + 1L])
}
R     <- get_arg("--reps", 50)
N_obs <- get_arg("--N", 400)
SEED  <- get_arg("--seed", 2026)

message(sprintf("Comparison study: R=%d, N=%d, seed=%d", R, N_obs, SEED))
set.seed(SEED)

# -- True parameters (in stanSPsurv parametrization) -------------------------

TRUE_BETA  <- c(-1.1, 0.45)
TRUE_GAMMA <- c(0.7, -0.5)
TRUE_SHAPE <- 0.75   # shape < 1 so BayesSPsurv can represent it

PARAM_NAMES <- c("beta[1]", "beta[2]", "gamma[1]", "gamma[2]", "shape")
TRUE_VALUES <- c(TRUE_BETA, TRUE_GAMMA, TRUE_SHAPE)

# -- Conditions to compare ---------------------------------------------------
# Pooled Weibull + Exchangeable Weibull (both within BayesSPsurv's support).
# Log-logistic omitted: BayesSPsurv's rho slice sampler for loglog uses a
# separate code path with the same case bug, and verification is pending.
# Spatial omitted: BayesSPsurv uses a different (unscaled) CAR prior.

conditions <- data.frame(
  structure = c("pooled", "exchangeable"),
  baseline  = c("weibull", "weibull"),
  stringsAsFactors = FALSE
)

# -- DGP ---------------------------------------------------------------------

simulate_comparison_data <- function(structure, baseline, N = N_obs, J = 8L) {
  x <- rnorm(N)
  z <- rnorm(N)
  X <- cbind(1, x)
  Z <- cbind(1, z)

  group_id <- rep(seq_len(J), length.out = N)

  # Frailty (survival stage only)
  u <- rep(0, J)
  if (structure == "exchangeable") {
    u <- rnorm(J, 0, 0.35)
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

  # Adjacency for spatial (not used in comparison conditions, but needed)
  A <- matrix(0, J, J)
  for (j in seq_len(J - 1L)) { A[j, j + 1L] <- 1; A[j + 1L, j] <- 1 }
  edges <- which(upper.tri(A) & A == 1, arr.ind = TRUE)
  graph <- stanSPsurv:::adjacency_graph(
    if (structure %in% c("spatial", "exchangeable")) A else NULL
  )

  # stanSPsurv data list
  stan_data <- list(
    N = N, M = N, K = ncol(X), Q = ncol(Z),
    t0 = rep(0, N), time = time, event = as.integer(event),
    X = X, Z_spell = Z,
    spell_start = seq_len(N), spell_end = seq_len(N),
    spell_event = as.integer(event),
    J = if (structure == "pooled") 1L else J,
    spell_group_id = if (structure == "pooled") rep(1L, N) else group_id,
    N_edges = if (structure != "pooled") nrow(edges) else 0L,
    node1 = if (structure != "pooled") as.integer(edges[, 1]) else integer(),
    node2 = if (structure != "pooled") as.integer(edges[, 2]) else integer(),
    N_components = graph$N_components,
    component_id = graph$component_id,
    component_size = graph$component_size,
    icar_scale = graph$scale,
    baseline = match(baseline, c("weibull", "loglogistic"))
  )

  # BayesSPsurv data.frame
  # C=1 means event occurred (failure/susceptible); C=0 means censored
  bayes_df <- data.frame(
    duration = time,
    t.0 = rep(0, N),
    cured = as.integer(event),  # C=1 if failure
    lastyear = rep(1L, N),      # single-row spells
    x = x,
    z = z,
    group_id = group_id
  )

  list(stan_data = stan_data, bayes_df = bayes_df, event_rate = mean(event),
       adjacency = A)
}

# -- Convert BayesSPsurv posteriors to stanSPsurv parametrization ------------

convert_bayes_to_stan <- function(betas_draws, gammas_draws, rho_draws,
                                  baseline) {
  # shape_stan = exp(-rho_bayes)
  shape_draws <- exp(-rho_draws)

  # beta_stan = -shape * beta_bayes
  beta_draws <- -sweep(betas_draws, 1, shape_draws, `*`)

  # gamma_stan = -gamma_bayes
  gamma_draws <- -gammas_draws

  list(beta = beta_draws, gamma = gamma_draws, shape = shape_draws)
}

# For log-logistic, BayesSPsurv uses rho=1 fixed (no shape sampling)
# and the parametrization is different. Let me check...
# Actually BayesSPsurv's loglog has rho fixed at 1 in the shape parameter
# and the baseline is: h(t) = r*lambda^r*t^(r-1) / (1 + (lambda*t)^r)
# with r = exp(-rho). For loglog in BayesSPsurv, rho IS sampled via
# a different slice sampler. Let me verify.

# -- Compile Stan models once ------------------------------------------------

message("Compiling Stan models...")
stan_models <- list()
for (s in unique(conditions$structure)) {
  stan_models[[s]] <- cmdstan_model(
    file.path("inst", "stan", paste0(s, "_spsurv.stan")),
    force_recompile = FALSE, quiet = TRUE
  )
}
message("Compilation complete.")

# -- Sampling settings -------------------------------------------------------

# stanSPsurv: 4 chains × warmup + sampling
STAN_CHAINS     <- 4L
STAN_WARMUP     <- as.integer(get_arg("--stan-warmup", 500))
STAN_SAMPLING   <- as.integer(get_arg("--stan-sampling", 500))

# BayesSPsurv: N total, burn, thin → (N - burn) / thin effective draws
# Default targets ~2000 draws to match Stan
BAYES_N    <- as.integer(get_arg("--bayes-N", 12000))
BAYES_BURN <- as.integer(get_arg("--bayes-burn", 2000))
BAYES_THIN <- as.integer(get_arg("--bayes-thin", 5))
# Effective draws: (N - burn) / thin

# -- Helper: summarize draws -------------------------------------------------

summarize_draws <- function(draws_matrix, param_names) {
  # draws_matrix: rows = draws, cols = parameters
  data.frame(
    param = param_names,
    mean  = colMeans(draws_matrix),
    q05   = apply(draws_matrix, 2, quantile, 0.05),
    q95   = apply(draws_matrix, 2, quantile, 0.95),
    ess   = apply(draws_matrix, 2, function(x) coda::effectiveSize(x)),
    stringsAsFactors = FALSE
  )
}

# -- Run comparison ----------------------------------------------------------

out_dir <- file.path("validation", "comparison_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- list()
idx <- 0L

for (cond in seq_len(nrow(conditions))) {
  structure <- conditions$structure[cond]
  baseline  <- conditions$baseline[cond]
  message(sprintf("\n=== %s / %s ===", structure, baseline))

  bayes_form <- if (baseline == "weibull") "weibull" else "loglog"

  for (r in seq_len(R)) {
    if (r %% 10 == 0 || r == 1) message(sprintf("  replication %d/%d", r, R))

    sim <- simulate_comparison_data(structure, baseline)

    # ---- stanSPsurv fit ----------------------------------------------------
    t_stan <- proc.time()
    stan_fit <- tryCatch(
      stan_models[[structure]]$sample(
        data = sim$stan_data,
        seed = SEED * 1000 + idx + 1L,
        chains = STAN_CHAINS, parallel_chains = STAN_CHAINS,
        iter_warmup = STAN_WARMUP, iter_sampling = STAN_SAMPLING,
        refresh = 0, adapt_delta = 0.90, max_treedepth = 11,
        show_messages = FALSE, show_exceptions = FALSE
      ),
      error = function(e) NULL
    )
    time_stan <- (proc.time() - t_stan)[["elapsed"]]

    # ---- BayesSPsurv fit ---------------------------------------------------
    t_bayes <- proc.time()
    bayes_fit <- tryCatch({
      if (structure == "pooled") {
        pooledSPsurv(
          duration = duration ~ x,
          immune   = cured ~ z,
          Y0       = "t.0",
          LY       = "lastyear",
          data     = sim$bayes_df,
          N        = BAYES_N,
          burn     = BAYES_BURN,
          thin     = BAYES_THIN,
          w        = c(1, 1, 1),
          m        = 10,
          form     = bayes_form
        )
      } else if (structure == "exchangeable") {
        exchangeSPsurv(
          duration  = duration ~ x,
          immune    = cured ~ z,
          Y0        = "t.0",
          LY        = "lastyear",
          S         = "group_id",
          data      = sim$bayes_df,
          N         = BAYES_N,
          burn      = BAYES_BURN,
          thin      = BAYES_THIN,
          w         = c(1, 1, 1),
          m         = 10,
          form      = bayes_form,
          prop.varV = 1e-04,
          prop.varW = 1e-04
        )
      }
    }, error = function(e) NULL)
    time_bayes <- (proc.time() - t_bayes)[["elapsed"]]

    # ---- Extract and convert stanSPsurv results ----------------------------
    stan_row <- NULL
    if (!is.null(stan_fit)) {
      draws <- stan_fit$draws(variables = c("beta", "gamma", "shape"),
                              format = "matrix")
      diag <- stan_fit$diagnostic_summary(quiet = TRUE)
      sel <- match(PARAM_NAMES, colnames(draws))
      d <- draws[, sel, drop = FALSE]
      stan_row <- data.frame(
        param = PARAM_NAMES, truth = TRUE_VALUES,
        mean = colMeans(d),
        q05 = apply(d, 2, quantile, 0.05),
        q95 = apply(d, 2, quantile, 0.95),
        ess = apply(d, 2, posterior::ess_bulk),
        rhat = apply(d, 2, posterior::rhat),
        divergences = sum(diag$num_divergent),
        stringsAsFactors = FALSE
      )
    }

    # ---- Extract and convert BayesSPsurv results --------------------------
    bayes_row <- NULL
    if (!is.null(bayes_fit)) {
      betas_draws  <- bayes_fit$betas   # matrix: draws × K
      gammas_draws <- bayes_fit$gammas   # matrix: draws × Q
      rho_draws    <- bayes_fit$rho      # vector: draws

      conv <- convert_bayes_to_stan(betas_draws, gammas_draws, rho_draws,
                                    baseline)

      # Combine into a matrix: [beta1, beta2, gamma1, gamma2, shape]
      draws_mat <- cbind(conv$beta, conv$gamma, conv$shape)
      colnames(draws_mat) <- PARAM_NAMES

      bsm <- summarize_draws(draws_mat, PARAM_NAMES)
      bayes_row <- data.frame(
        param = PARAM_NAMES, truth = TRUE_VALUES,
        mean = bsm$mean, q05 = bsm$q05, q95 = bsm$q95,
        ess = bsm$ess, rhat = NA_real_,
        divergences = NA_integer_,
        stringsAsFactors = FALSE
      )
    }

    idx <- idx + 1L

    if (!is.null(stan_row)) {
      all_results[[length(all_results) + 1]] <- cbind(
        data.frame(structure = structure, baseline = baseline, rep = r,
                   package = "stanSPsurv", elapsed_sec = time_stan,
                   event_rate = sim$event_rate, stringsAsFactors = FALSE),
        stan_row
      )
    }
    if (!is.null(bayes_row)) {
      all_results[[length(all_results) + 1]] <- cbind(
        data.frame(structure = structure, baseline = baseline, rep = r,
                   package = "BayesSPsurv", elapsed_sec = time_bayes,
                   event_rate = sim$event_rate, stringsAsFactors = FALSE),
        bayes_row
      )
    }
  }
}

results <- do.call(rbind, all_results)
write.csv(results, file.path(out_dir, "raw_comparison.csv"), row.names = FALSE)

# -- Summary table -----------------------------------------------------------

message("\nComputing summary statistics...")

summary_table <- do.call(rbind, lapply(
  split(results,
        interaction(results$structure, results$baseline,
                    results$package, results$param)),
  function(df) {
    df <- df[!is.na(df$mean), ]
    if (nrow(df) == 0) return(NULL)
    bias     <- mean(df$mean - df$truth)
    rmse     <- sqrt(mean((df$mean - df$truth)^2))
    coverage <- mean(df$q05 <= df$truth & df$truth <= df$q95)
    med_ess  <- median(df$ess, na.rm = TRUE)
    med_time <- median(df$elapsed_sec)
    ess_sec  <- med_ess / med_time

    data.frame(
      structure = df$structure[1], baseline = df$baseline[1],
      package = df$package[1], param = df$param[1], truth = df$truth[1],
      bias = bias, rmse = rmse, coverage_90 = coverage,
      median_ess = med_ess, ess_per_sec = ess_sec,
      median_seconds = med_time, n_reps = nrow(df),
      stringsAsFactors = FALSE
    )
  }
))
rownames(summary_table) <- NULL
write.csv(summary_table, file.path(out_dir, "summary_comparison.csv"),
          row.names = FALSE)

message("\n=== COMPARISON SUMMARY ===")
print(summary_table, digits = 3, row.names = FALSE)
message(sprintf("\nResults saved to %s/", out_dir))
