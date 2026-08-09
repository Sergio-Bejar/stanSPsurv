/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
#' Fit a pooled split-population survival model
#'
#' @param formula Survival-stage formula with a two-column response
#'   `Surv(time, event)`.
#' @param split_formula Formula for the latent susceptibility probability.
#' @param data A data frame.
#' @param spell Episode identifier. Required for counting-process responses;
#'   split-stage covariates must be constant within each episode.
#' @param form Event-time distribution: `"weibull"` or `"loglogistic"`.
#' @param ... Sampling arguments passed to [cmdstanr::CmdStanModel$sample()].
#' @return An object of class `stanSPsurv`.
#' @export
stan_pooled_spsurv <- function(formula, split_formula, data, spell = NULL,
                               form = c("weibull", "loglogistic"), ...) {
  spell <- resolve_optional(substitute(spell), spell, data, parent.frame())
  fit_spsurv(formula, split_formula, data, form, structure = "pooled",
             spell = spell, ...)
}

#' Fit an exchangeable-frailty split-population survival model
#' @inheritParams stan_pooled_spsurv
#' @param group Group identifier, supplied as a column name or vector.
#' @export
stan_exchange_spsurv <- function(formula, split_formula, data, group,
                                 spell = NULL,
                                 form = c("weibull", "loglogistic"), ...) {
  group <- resolve_group(substitute(group), data, parent.frame())
  spell <- resolve_optional(substitute(spell), spell, data, parent.frame())
  fit_spsurv(formula, split_formula, data, form, structure = "exchangeable",
             group = group, spell = spell, ...)
}

#' Fit a spatial CAR split-population survival model
#' @inheritParams stan_pooled_spsurv
#' @param area Spatial unit identifier, supplied as a column name or vector.
#' @param adjacency Symmetric binary adjacency matrix ordered by the unique
#'   values of `area`.
#' @export
stan_spatial_spsurv <- function(formula, split_formula, data, area, adjacency,
                                spell = NULL,
                                form = c("weibull", "loglogistic"), ...) {
  area <- resolve_group(substitute(area), data, parent.frame())
  spell <- resolve_optional(substitute(spell), spell, data, parent.frame())
  validate_adjacency(adjacency, length(unique(area)))
  fit_spsurv(formula, split_formula, data, form, structure = "spatial",
             group = area, adjacency = adjacency, spell = spell, ...)
}

fit_spsurv <- function(formula, split_formula, data, form, structure,
                       group = NULL, adjacency = NULL, spell = NULL, ...) {
  form <- match.arg(tolower(form), c("weibull", "loglogistic"))
  prep <- prepare_spsurv_data(formula, split_formula, data, group, adjacency,
                              spell = spell)
  prep$stan_data$baseline <- match(form, c("weibull", "loglogistic"))
  stan_name <- paste0(structure, "_spsurv.stan")
  stan_file <- system.file("stan", stan_name, package = "stanSPsurv")
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

resolve_optional <- function(expr, value, data, env) {
  if (identical(expr, quote(NULL)) || is.null(value)) return(NULL)
  resolve_group(expr, data, env)
}

resolve_group <- function(expr, data, env) {
  value <- try(eval(expr, data, env), silent = TRUE)
  if (inherits(value, "try-error")) stop("Could not evaluate the group/area variable.", call. = FALSE)
  if (length(value) != nrow(data)) stop("Group/area must have one value per row.", call. = FALSE)
  value
}
