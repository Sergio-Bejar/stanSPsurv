/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
#!/usr/bin/env Rscript

# Short smoke/recovery study for all six candidate model specifications.
# These runs establish compilation and basic sampling behavior; they do not
# replace a larger simulation study or comparison with BayesSPsurv.

local_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

cmdstan_dir <- file.path(getwd(), ".cmdstan", "cmdstan-2.39.0")
if (dir.exists(cmdstan_dir)) cmdstanr::set_cmdstan_path(cmdstan_dir)

set.seed(20260808)
dir.create(file.path("validation", "output"), recursive = TRUE, showWarnings = FALSE)

simulate_case <- function(structure, baseline, N = 180L, J = 6L) {
  x <- rnorm(N)
  z <- rnorm(N)
  X <- cbind(1, x)
  Z <- cbind(1, z)
  beta <- c(-1.1, 0.45)
  gamma <- c(0.7, -0.5)
  shape <- 1.35
  group_id <- rep(seq_len(J), length.out = N)
  u <- rep(0, J)
  if (structure == "exchangeable") {
    u <- seq(-0.45, 0.45, length.out = J)
  } else if (structure == "spatial") {
    u <- c(-0.45, -0.25, -0.05, 0.10, 0.25, 0.40)
    u <- u - mean(u)
  }
  pi <- plogis(drop(Z %*% gamma))
  susceptible <- rbinom(N, 1, pi)
  eta <- drop(X %*% beta) + u[group_id]
  q <- runif(N)
  if (baseline == "weibull") {
    latent <- (-log1p(-q) / exp(eta))^(1 / shape)
  } else {
    latent <- (q / ((1 - q) * exp(eta)))^(1 / shape)
  }
  latent[susceptible == 0] <- Inf
  censor <- rexp(N, rate = 0.10)
  time <- pmin(latent, censor)
  event <- as.integer(latent <= censor)
  A <- matrix(0, J, J)
  for (j in seq_len(J - 1L)) A[j, j + 1L] <- A[j + 1L, j] <- 1
  edges <- which(upper.tri(A) & A == 1, arr.ind = TRUE)
  graph <- stanSPsurv:::adjacency_graph(if (structure == "spatial") A else NULL)
  list(
    data = list(N = N, M = N, K = ncol(X), Q = ncol(Z), time = time,
      t0 = rep(0, N), event = event, X = X, Z_spell = Z,
      spell_start = seq_len(N), spell_end = seq_len(N), spell_event = event,
      J = if (structure == "pooled") 1L else J,
      spell_group_id = if (structure == "pooled") rep(1L, N) else group_id,
      N_edges = if (structure == "spatial") nrow(edges) else 0L,
      node1 = if (structure == "spatial") as.integer(edges[, 1]) else integer(),
      node2 = if (structure == "spatial") as.integer(edges[, 2]) else integer(),
      N_components = graph$N_components, component_id = graph$component_id,
      component_size = graph$component_size, icar_scale = graph$scale,
      baseline = match(baseline, c("weibull", "loglogistic"))),
    truth = c(beta_1 = beta[1], beta_2 = beta[2], gamma_1 = gamma[1],
              gamma_2 = gamma[2], shape = shape),
    event_rate = mean(event)
  )
}

models <- lapply(c("pooled", "exchangeable", "spatial"), function(s)
  cmdstanr::cmdstan_model(file.path("inst", "stan", paste0(s, "_spsurv.stan")),
                         force_recompile = TRUE, quiet = TRUE))
names(models) <- c("pooled", "exchangeable", "spatial")
cases <- expand.grid(
  structure = c("pooled", "exchangeable", "spatial"),
  baseline = c("weibull", "loglogistic"),
  stringsAsFactors = FALSE
)

results <- vector("list", nrow(cases))
for (i in seq_len(nrow(cases))) {
  structure <- cases$structure[i]
  baseline <- cases$baseline[i]
  message("Running ", structure, " / ", baseline)
  sim <- simulate_case(structure, baseline)
  case_dir <- file.path("validation", "output", paste(structure, baseline, sep = "-"))
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- models[[structure]]$sample(
    data = sim$data, seed = 1000 + i, chains = 2, parallel_chains = 2,
    iter_warmup = 300, iter_sampling = 300, refresh = 0,
    adapt_delta = 0.90, max_treedepth = 11,
    output_dir = case_dir
  )
  variables <- c("beta", "gamma", "shape")
  if (structure != "pooled") variables <- c(variables, "sigma_v", "sigma_w")
  sm <- fit$summary(variables = variables)
  diag <- fit$diagnostic_summary()
  wanted <- c("beta[1]", "beta[2]", "gamma[1]", "gamma[2]", "shape")
  est <- sm$mean[match(wanted, sm$variable)]
  results[[i]] <- data.frame(
    structure = structure, baseline = baseline, event_rate = sim$event_rate,
    max_rhat = max(sm$rhat, na.rm = TRUE), min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
    divergences = sum(diag$num_divergent), max_treedepth = sum(diag$num_max_treedepth),
    beta_1_error = est[1] - sim$truth[1], beta_2_error = est[2] - sim$truth[2],
    gamma_1_error = est[3] - sim$truth[3], gamma_2_error = est[4] - sim$truth[4],
    shape_error = est[5] - sim$truth[5]
  )
}

results <- do.call(rbind, results)
write.csv(results, file.path("validation", "results.csv"), row.names = FALSE)
print(results, digits = 3, row.names = FALSE)
