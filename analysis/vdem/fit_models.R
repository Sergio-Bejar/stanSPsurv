#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
structure <- if (length(args) >= 1) tolower(args[[1]]) else "pooled"
form <- if (length(args) >= 2) tolower(args[[2]]) else "weibull"
if (!structure %in% c("pooled", "exchangeable", "spatial"))
  stop("First argument must be pooled, exchangeable, or spatial.")
if (!form %in% c("weibull", "loglogistic"))
  stop("Second argument must be weibull or loglogistic.")

if (!requireNamespace("survival", quietly = TRUE) ||
    !requireNamespace("stanSPsurv", quietly = TRUE) ||
    !requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Install survival, cmdstanr, and this stanSPsurv package before fitting.")
}

root <- normalizePath(getwd(), mustWork = TRUE)
input <- file.path(root, "analysis", "vdem", "output", "vdem_application.rds")
if (!file.exists(input)) stop("Run analysis/vdem/prepare_vdem.R first.")
app <- readRDS(input)
d <- app$data

survival_formula <- survival::Surv(duration, breakdown) ~
  judicial_constraints_z + legislative_constraints_z + free_expression_z +
  log_gdppc_z + civil_war_z + prior_dem_years_z
split_formula <- ~ log_gdppc_z + prior_dem_years_z

sampling <- list(chains = 4, parallel_chains = 4, iter_warmup = 1000,
                 iter_sampling = 1000, seed = 16032026,
                 adapt_delta = 0.97, max_treedepth = 13,
                 refresh = 100)
if (structure == "pooled") {
  fit <- do.call(stanSPsurv::stan_pooled_spsurv,
                 c(list(survival_formula, split_formula, d, form = form), sampling))
} else if (structure == "exchangeable") {
  fit <- do.call(stanSPsurv::stan_exchange_spsurv,
                 c(list(survival_formula, split_formula, d,
                        group = d$country_text_id, form = form), sampling))
} else {
  fit <- do.call(stanSPsurv::stan_spatial_spsurv,
                 c(list(survival_formula, split_formula, d,
                        area = d$country_text_id, adjacency = app$adjacency,
                        form = form), sampling))
}

fit_dir <- file.path(root, "analysis", "vdem", "output", "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path(fit_dir, paste0(structure, "_", form, ".rds")))
model_summary <- summary(fit)
utils::write.csv(model_summary,
                 file.path(root, "analysis", "vdem", "output",
                           paste0(structure, "_", form, "_summary.csv")),
                 row.names = FALSE)
diagnostics <- fit$fit$diagnostic_summary()
diagnostic_table <- data.frame(
  chain = seq_along(diagnostics$num_divergent),
  divergent = diagnostics$num_divergent,
  max_treedepth = diagnostics$num_max_treedepth,
  ebfmi = diagnostics$ebfmi
)
utils::write.csv(diagnostic_table,
                 file.path(root, "analysis", "vdem", "output",
                           paste0(structure, "_", form, "_diagnostics.csv")),
                 row.names = FALSE)
print(model_summary)
