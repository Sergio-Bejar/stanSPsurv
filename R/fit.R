#' Fit a pooled split-population survival model
#'
#' @param formula Survival-stage formula with a two-column response
#'   `Surv(time, event)`.
#' @param split_formula Formula for the latent susceptibility probability.
#' @param data A data frame.
#' @param form Event-time distribution: `"weibull"` or `"loglogistic"`.
#' @param ... Sampling arguments passed to [cmdstanr::CmdStanModel$sample()].
#' @return An object of class `stanSPsurv`.
#' @export
stan_pooled_spsurv <- function(formula, split_formula, data,
                               form = c("weibull", "loglogistic"), ...) {
  fit_spsurv(formula, split_formula, data, form, structure = "pooled", ...)
}

#' Fit an exchangeable-frailty split-population survival model
#' @inheritParams stan_pooled_spsurv
#' @param group Group identifier, supplied as a column name or vector.
#' @export
stan_exchange_spsurv <- function(formula, split_formula, data, group,
                                 form = c("weibull", "loglogistic"), ...) {
  group <- resolve_group(substitute(group), data, parent.frame())
  fit_spsurv(formula, split_formula, data, form, structure = "exchangeable",
             group = group, ...)
}

#' Fit a spatial CAR split-population survival model
#' @inheritParams stan_pooled_spsurv
#' @param area Spatial unit identifier, supplied as a column name or vector.
#' @param adjacency Symmetric binary adjacency matrix ordered by the unique
#'   values of `area`.
#' @export
stan_spatial_spsurv <- function(formula, split_formula, data, area, adjacency,
                                form = c("weibull", "loglogistic"), ...) {
  area <- resolve_group(substitute(area), data, parent.frame())
  validate_adjacency(adjacency, length(unique(area)))
  fit_spsurv(formula, split_formula, data, form, structure = "spatial",
             group = area, adjacency = adjacency, ...)
}

fit_spsurv <- function(formula, split_formula, data, form, structure,
                       group = NULL, adjacency = NULL, ...) {
  form <- match.arg(tolower(form), c("weibull", "loglogistic"))
  prep <- prepare_spsurv_data(formula, split_formula, data, group, adjacency)
  prep$stan_data$baseline <- match(form, c("weibull", "loglogistic"))
  prep$stan_data$structure_type <- match(structure,
                                        c("pooled", "exchangeable", "spatial")) - 1L
  stan_file <- system.file("stan", "spsurv.stan", package = "stanSPsurv")
  if (!nzchar(stan_file)) stop("Installed Stan model was not found.", call. = FALSE)
  model <- cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
  fit <- model$sample(data = prep$stan_data, ...)
  out <- list(fit = fit, formula = formula, split_formula = split_formula,
              form = form, structure = structure, terms = prep$terms,
              split_terms = prep$split_terms, group_levels = prep$group_levels,
              call = match.call())
  class(out) <- "stanSPsurv"
  out
}

resolve_group <- function(expr, data, env) {
  value <- try(eval(expr, data, env), silent = TRUE)
  if (inherits(value, "try-error")) stop("Could not evaluate the group/area variable.", call. = FALSE)
  if (length(value) != nrow(data)) stop("Group/area must have one value per row.", call. = FALSE)
  value
}

