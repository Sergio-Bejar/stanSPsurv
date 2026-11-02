#' @export
print.stanSPsurv <- function(x, ...) {
  cat("stanSPsurv fit\n", "  structure: ", x$structure,
      "\n  baseline:  ", x$form, "\n", sep = "")
  invisible(x)
}

#' @export
summary.stanSPsurv <- function(object, ...) {
  variables <- c("beta", "gamma", "shape")
  if (object$structure != "pooled") variables <- c(variables, "sigma_v", "sigma_w")
  if (object$structure == "spatial") variables <- c(variables, "phi_v", "phi_w")
  object$fit$summary(variables = variables)
}

#' @export
coef.stanSPsurv <- function(object, component = c("survival", "split"), ...) {
  component <- match.arg(component)
  variable <- if (component == "survival") "beta" else "gamma"
  object$fit$summary(variables = variable)
}

#' Extract group or spatial frailties
#' @param object A fitted `stanSPsurv` object.
#' @param ... Reserved for extensions.
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

#' @export
predict.stanSPsurv <- function(object, newdata = NULL,
                               type = c("susceptibility", "survival"),
                               times = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) stop("newdata is required in this development version.", call. = FALSE)
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
