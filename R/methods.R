#' Print a stanSPsurv fit
#'
#' @param x A fitted `stanSPsurv` object.
#' @param ... Not used.
#' @return Invisible `x`.
#' @export
print.stanSPsurv <- function(x, ...) {
  cat("stanSPsurv fit\n", "  structure: ", x$structure,
      "\n  baseline:  ", x$form, "\n", sep = "")
  invisible(x)
}

#' Summarise a stanSPsurv fit
#'
#' Returns posterior summaries (mean, median, sd, quantiles, Rhat, ESS)
#' for all model parameters.
#'
#' @param object A fitted `stanSPsurv` object.
#' @param ... Not used.
#' @return A tibble of posterior summaries produced by
#'   [cmdstanr::CmdStanMCMC].
#' @export
summary.stanSPsurv <- function(object, ...) {
  variables <- c("beta", "gamma", "shape")
  if (object$structure != "pooled") variables <- c(variables, "sigma_v", "sigma_w")
  if (object$structure == "spatial") variables <- c(variables, "phi_v", "phi_w")
  object$fit$summary(variables = variables)
}

#' Extract coefficients from a stanSPsurv fit
#'
#' @param object A fitted `stanSPsurv` object.
#' @param component Which stage: `"survival"` (beta) or `"split"` (gamma).
#' @param ... Not used.
#' @return A tibble of posterior summaries for the requested coefficients.
#' @export
coef.stanSPsurv <- function(object, component = c("survival", "split"), ...) {
  component <- match.arg(component)
  variable <- if (component == "survival") "beta" else "gamma"
  object$fit$summary(variables = variable)
}

#' Extract group or spatial frailties
#'
#' Returns posterior summaries for the group-level random effects
#' (frailties) from an exchangeable or spatial model.
#'
#' @param object A fitted `stanSPsurv` object.
#' @param stage Which stage: `"survival"` (w) or `"split"` (v).
#' @param ... Not used.
#' @return A data frame with posterior summaries and group labels.
#'   Returns an empty data frame for pooled models.
#' @export
frailty <- function(object, stage = c("survival", "split"), ...) {
  if (!inherits(object, "stanSPsurv")) stop("object must be a stanSPsurv fit.", call. = FALSE)
  if (object$structure == "pooled") return(data.frame())
  stage <- match.arg(stage)
  variable <- if (stage == "survival") "w" else "v"
  ans <- object$fit$summary(variables = variable)
  ans$group <- object$group_levels
  ans$stage <- stage
  ans
}

#' Predict from a stanSPsurv fit
#'
#' Compute posterior-mean susceptibility probabilities or survival
#' probabilities for new covariate profiles.
#'
#' @param object A fitted `stanSPsurv` object.
#' @param newdata A data frame of covariate profiles.
#' @param type `"susceptibility"` for P(susceptible) or `"survival"`
#'   for S(t).
#' @param times Numeric vector of positive time points (required when
#'   `type = "survival"`).
#' @param ... Not used.
#' @return For `"susceptibility"`, a numeric vector of posterior-mean
#'   probabilities (one per row of `newdata`).
#'   For `"survival"`, a matrix of posterior-mean survival probabilities
#'   (rows = profiles, columns = times).
#' @export
predict.stanSPsurv <- function(object, newdata = NULL,
                               type = c("susceptibility", "survival"),
                               times = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) stop("newdata is required.", call. = FALSE)
  draws <- object$fit$draws(c("beta", "gamma", "shape"), format = "matrix")
  if (type == "susceptibility") {
    Z <- stats::model.matrix(stats::delete.response(object$split_terms), newdata)
    g <- draws[, grep("^gamma\\[", colnames(draws)), drop = FALSE]
    return(colMeans(stats::plogis(g %*% t(Z))))
  }
  if (is.null(times) || any(times <= 0)) stop("Positive times are required.", call. = FALSE)
  X <- stats::model.matrix(stats::delete.response(object$terms), newdata)
  b <- draws[, grep("^beta\\[", colnames(draws)), drop = FALSE]
  shape <- draws[, "shape"]
  eta <- b %*% t(X)
  out <- vapply(times, function(tt) {
    S <- if (object$form == "weibull") exp(-exp(eta) * tt^shape) else 1 / (1 + exp(eta) * tt^shape)
    colMeans(S)
  }, numeric(nrow(X)))
  if (length(times) == 1L) drop(out) else out
}
