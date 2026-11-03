#' Plot a stanSPsurv fit
#'
#' Produces a coefficient forest plot showing posterior means and 90\%
#' credible intervals for the survival-stage (`beta`) and split-stage
#' (`gamma`) coefficients. Intercepts are omitted by default.
#'
#' @param x A fitted `stanSPsurv` object.
#' @param which Which components to plot: `"both"`, `"survival"`, or
#'   `"split"`.
#' @param include_intercept Logical; include intercept terms?
#' @param prob Width of the credible interval (default 0.90).
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot.default()].
#' @return Invisible `NULL`. Called for its side effect.
#' @export
plot.stanSPsurv <- function(x, which = c("both", "survival", "split"),
                            include_intercept = FALSE, prob = 0.90,
                            main = "Coefficient plot", ...) {
  which <- match.arg(which)
  alpha <- (1 - prob) / 2
  draws <- x$fit$draws(c("beta", "gamma"), format = "matrix")

  beta_idx  <- grep("^beta\\[", colnames(draws))
  gamma_idx <- grep("^gamma\\[", colnames(draws))

  beta_labs  <- coef_labels(x$terms)
  gamma_labs <- coef_labels(x$split_terms)

  sections <- list()
  if (which %in% c("both", "survival"))
    sections$survival <- list(draws = draws[, beta_idx, drop = FALSE],
                              labs = beta_labs, title = "Survival (beta)")
  if (which %in% c("both", "split"))
    sections$split <- list(draws = draws[, gamma_idx, drop = FALSE],
                           labs = gamma_labs, title = "Split (gamma)")

  # Drop intercepts unless requested
  if (!include_intercept) {
    for (nm in names(sections)) {
      keep <- !grepl("^\\(Intercept\\)$", sections[[nm]]$labs)
      sections[[nm]]$draws <- sections[[nm]]$draws[, keep, drop = FALSE]
      sections[[nm]]$labs  <- sections[[nm]]$labs[keep]
    }
  }

  all_labs <- unlist(lapply(rev(sections), `[[`, "labs"))
  n <- length(all_labs)
  if (n == 0L) { message("No coefficients to plot."); return(invisible(NULL)) }

  means <- numeric(n); lo <- numeric(n); hi <- numeric(n)
  grp <- character(n)
  pos <- 0L
  for (sec in rev(sections)) {
    for (j in seq_along(sec$labs)) {
      pos <- pos + 1L
      d <- sec$draws[, j]
      means[pos] <- mean(d)
      lo[pos]    <- stats::quantile(d, alpha)
      hi[pos]    <- stats::quantile(d, 1 - alpha)
      grp[pos]   <- sec$title
    }
  }

  xlim <- range(c(lo, hi, 0))
  xlim <- xlim + diff(xlim) * c(-0.1, 0.1)
  cols <- ifelse(grp == "Survival (beta)", "steelblue", "firebrick")
  pch  <- ifelse(grp == "Survival (beta)", 19, 17)

  opar <- graphics::par(mar = c(4, max(nchar(all_labs)) * 0.55 + 1, 2, 1))
  on.exit(graphics::par(opar))

  plot(means, seq_len(n), xlim = xlim, yaxt = "n", pch = pch, col = cols,
       xlab = "Estimate", ylab = "", main = main,
       panel.first = graphics::abline(v = 0, lty = 2, col = "grey60"), ...)
  graphics::segments(lo, seq_len(n), hi, seq_len(n), col = cols, lwd = 2)
  graphics::axis(2, at = seq_len(n), labels = all_labs, las = 1, cex.axis = 0.8)

  ugrp <- unique(grp)
  if (length(ugrp) > 1L)
    graphics::legend("topright", legend = ugrp,
                     col = c("steelblue", "firebrick")[seq_along(ugrp)],
                     pch = c(19, 17)[seq_along(ugrp)], lwd = 2, bty = "n",
                     cex = 0.8)
  invisible(NULL)
}

#' Plot survival curves from a stanSPsurv fit
#'
#' Draws posterior-mean survival curves (and pointwise credible bands)
#' for one or more covariate profiles.
#'
#' @param object A fitted `stanSPsurv` object.
#' @param newdata A data frame of covariate profiles (one row per curve).
#' @param times Numeric vector of time points at which to evaluate survival.
#' @param prob Width of the credible band (default 0.90).
#' @param labels Optional character vector of curve labels (length
#'   `nrow(newdata)`).
#' @param cols Optional colour vector.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot.default()].
#' @return Invisible matrix of posterior-mean survival values
#'   (rows = profiles, cols = times).
#' @export
plot_survival <- function(object, newdata, times, prob = 0.90,
                          labels = NULL, cols = NULL,
                          main = "Survival curves", ...) {
  if (!inherits(object, "stanSPsurv"))
    stop("object must be a stanSPsurv fit.", call. = FALSE)
  if (missing(times) || length(times) < 2L)
    stop("Provide a numeric vector of at least two time points.", call. = FALSE)

  alpha <- (1 - prob) / 2
  draws <- object$fit$draws(c("beta", "shape"), format = "matrix")
  b     <- draws[, grep("^beta\\[", colnames(draws)), drop = FALSE]
  shape <- draws[, "shape"]
  X     <- stats::model.matrix(stats::delete.response(object$terms), newdata)
  np    <- nrow(X)
  nt    <- length(times)

  S_mean <- matrix(NA_real_, np, nt)
  S_lo   <- matrix(NA_real_, np, nt)
  S_hi   <- matrix(NA_real_, np, nt)

  for (i in seq_len(np)) {
    eta_i <- as.vector(b %*% X[i, ])
    for (j in seq_len(nt)) {
      tt <- times[j]
      Sij <- if (object$form == "weibull") {
        exp(-exp(eta_i) * tt^shape)
      } else {
        1 / (1 + exp(eta_i) * tt^shape)
      }
      S_mean[i, j] <- mean(Sij)
      S_lo[i, j]   <- stats::quantile(Sij, alpha)
      S_hi[i, j]   <- stats::quantile(Sij, 1 - alpha)
    }
  }

  if (is.null(cols))
    cols <- if (np <= 8) seq_len(np) + 1L else grDevices::rainbow(np)
  if (is.null(labels)) labels <- paste0("Profile ", seq_len(np))

  plot(range(times), c(0, 1), type = "n", xlab = "Time",
       ylab = "S(t)", main = main, ...)
  for (i in seq_len(np)) {
    graphics::polygon(c(times, rev(times)),
                      c(S_lo[i, ], rev(S_hi[i, ])),
                      col = grDevices::adjustcolor(cols[i], alpha.f = 0.15),
                      border = NA)
    graphics::lines(times, S_mean[i, ], col = cols[i], lwd = 2)
  }
  graphics::legend("topright", legend = labels, col = cols, lwd = 2,
                   bty = "n", cex = 0.8)
  invisible(S_mean)
}

#' Plot susceptibility across a covariate range
#'
#' Shows how the posterior-mean probability of susceptibility varies
#' with one split-stage covariate, holding others at their means.
#'
#' @param object A fitted `stanSPsurv` object.
#' @param variable Name of the split-formula covariate to vary.
#' @param data The original data frame used for fitting (needed
#'   to compute covariate means for held-constant variables).
#' @param n_points Number of points along the covariate range (default 100).
#' @param prob Width of the credible band (default 0.90).
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot.default()].
#' @return Invisible data frame with columns `x`, `mean`, `lo`, `hi`.
#' @export
plot_susceptibility <- function(object, variable, data, n_points = 100L,
                                prob = 0.90,
                                main = "Susceptibility", ...) {
  if (!inherits(object, "stanSPsurv"))
    stop("object must be a stanSPsurv fit.", call. = FALSE)
  if (!variable %in% names(data))
    stop(paste0("'", variable, "' not found in data."), call. = FALSE)

  alpha <- (1 - prob) / 2
  x_range <- range(data[[variable]], na.rm = TRUE)
  x_seq   <- seq(x_range[1], x_range[2], length.out = n_points)

  # Build newdata: vary target, hold others at mean
  newdata <- as.data.frame(lapply(data, function(col) {
    if (is.numeric(col)) rep(mean(col, na.rm = TRUE), n_points) else rep(col[1], n_points)
  }))
  newdata[[variable]] <- x_seq

  draws <- object$fit$draws("gamma", format = "matrix")
  g     <- draws[, grep("^gamma\\[", colnames(draws)), drop = FALSE]
  Z     <- stats::model.matrix(stats::delete.response(object$split_terms), newdata)
  logit_pi <- g %*% t(Z)   # draws × n_points

  pi_mean <- colMeans(stats::plogis(logit_pi))
  pi_lo   <- apply(stats::plogis(logit_pi), 2, stats::quantile, alpha)
  pi_hi   <- apply(stats::plogis(logit_pi), 2, stats::quantile, 1 - alpha)

  plot(x_seq, pi_mean, type = "l", lwd = 2, col = "steelblue",
       xlab = variable, ylab = expression(P(susceptible)),
       ylim = c(0, 1), main = main, ...)
  graphics::polygon(c(x_seq, rev(x_seq)), c(pi_lo, rev(pi_hi)),
                    col = grDevices::adjustcolor("steelblue", alpha.f = 0.2),
                    border = NA)
  graphics::lines(x_seq, pi_mean, lwd = 2, col = "steelblue")

  invisible(data.frame(x = x_seq, mean = pi_mean, lo = pi_lo, hi = pi_hi))
}

#' Plot group or spatial frailties
#'
#' Caterpillar plot of posterior frailty estimates with credible intervals,
#' ordered by posterior mean.
#'
#' @param object A fitted `stanSPsurv` object (exchangeable or spatial).
#' @param stage Which frailty: `"survival"` (w) or `"split"` (v).
#' @param prob Width of the credible interval (default 0.90).
#' @param ... Additional arguments passed to [plot.default()].
#' @return Invisible data frame of frailty summaries.
#' @export
plot_frailty <- function(object, stage = c("survival", "split"),
                         prob = 0.90, ...) {
  if (!inherits(object, "stanSPsurv"))
    stop("object must be a stanSPsurv fit.", call. = FALSE)
  if (object$structure == "pooled")
    stop("Frailty plots require an exchangeable or spatial model.", call. = FALSE)

  stage <- match.arg(stage)
  alpha <- (1 - prob) / 2
  variable <- if (stage == "survival") "w" else "v"
  draws <- object$fit$draws(variable, format = "matrix")
  J <- ncol(draws)
  labs <- if (!is.null(object$group_levels)) object$group_levels else seq_len(J)

  means <- colMeans(draws)
  lo    <- apply(draws, 2, stats::quantile, alpha)
  hi    <- apply(draws, 2, stats::quantile, 1 - alpha)

  ord <- order(means)
  means <- means[ord]; lo <- lo[ord]; hi <- hi[ord]; labs <- labs[ord]

  stage_label <- if (stage == "survival") "Survival (w)" else "Split (v)"

  opar <- graphics::par(mar = c(4, max(nchar(as.character(labs))) * 0.55 + 1, 2, 1))
  on.exit(graphics::par(opar))

  xlim <- range(c(lo, hi, 0))
  xlim <- xlim + diff(xlim) * c(-0.1, 0.1)
  cols <- ifelse(lo > 0, "steelblue", ifelse(hi < 0, "firebrick", "grey40"))

  plot(means, seq_len(J), xlim = xlim, yaxt = "n", pch = 19, col = cols,
       xlab = "Frailty", ylab = "",
       main = paste("Frailties:", stage_label),
       panel.first = graphics::abline(v = 0, lty = 2, col = "grey60"), ...)
  graphics::segments(lo, seq_len(J), hi, seq_len(J), col = cols, lwd = 2)
  graphics::axis(2, at = seq_len(J), labels = labs, las = 1, cex.axis = 0.7)
  invisible(data.frame(group = labs, mean = means, lo = lo, hi = hi))
}

#' Choropleth map of spatial frailties
#'
#' Draws a choropleth (filled map) of posterior-mean frailties from a
#' spatial `stanSPsurv` model.  Two back-ends are supported:
#'
#' * **`maps` back-end** (default): pass `regions`, a character vector of
#'   country or region names matching the `maps::world` database, in the
#'   same order as the spatial units used when fitting.
#' * **`sf` back-end**: pass `sf_map`, an `sf` object with a geometry
#'   column and a key column whose values match `group_levels`.
#'
#' @param object A fitted `stanSPsurv` object with `structure = "spatial"`.
#' @param stage Which frailty: `"survival"` (w) or `"split"` (v).
#' @param regions Character vector of region names for the `maps` back-end,
#'   ordered to match `group_levels`.
#' @param sf_map An `sf` object for the `sf` back-end.
#' @param sf_key Column name in `sf_map` that matches `group_levels`.
#' @param palette A vector of colours for the gradient.  Default is a
#'   blue-white-red diverging palette.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot()].
#' @return Invisible data frame with group, mean frailty, and assigned colour.
#' @export
plot_choropleth <- function(object,
                            stage = c("survival", "split"),
                            regions = NULL,
                            sf_map = NULL, sf_key = NULL,
                            palette = NULL,
                            main = NULL, ...) {
  if (!inherits(object, "stanSPsurv"))
    stop("object must be a stanSPsurv fit.", call. = FALSE)
  if (object$structure != "spatial")
    stop("Choropleth requires a spatial model.", call. = FALSE)

  stage <- match.arg(stage)
  variable <- if (stage == "survival") "w" else "v"
  draws <- object$fit$draws(variable, format = "matrix")
  means <- colMeans(draws)
  labs  <- object$group_levels
  if (is.null(labs)) labs <- seq_along(means)

  if (is.null(main))
    main <- paste0("Spatial frailties: ", if (stage == "survival") "w" else "v")

  # Colour ramp
  if (is.null(palette))
    palette <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
                 "#FDDBC7", "#EF8A62", "#B2182B")
  ramp <- grDevices::colorRampPalette(palette)
  n_cols <- 100
  col_seq <- ramp(n_cols)

  # Map values to colours (diverging around 0)
  abs_max <- max(abs(means))
  if (abs_max == 0) abs_max <- 1
  scaled <- (means / abs_max + 1) / 2   # rescale to [0, 1]
  col_idx <- pmax(1L, pmin(n_cols, as.integer(scaled * (n_cols - 1)) + 1L))
  fill_cols <- col_seq[col_idx]

  # ---- sf back-end ---------------------------------------------------------
  if (!is.null(sf_map)) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required for the sf back-end.", call. = FALSE)
    if (is.null(sf_key))
      stop("sf_key must name the column in sf_map that matches group_levels.",
           call. = FALSE)
    sf_map$.frailty_mean <- means[match(sf_map[[sf_key]], labs)]
    sf_map$.frailty_col  <- fill_cols[match(sf_map[[sf_key]], labs)]
    plot(sf::st_geometry(sf_map), col = sf_map$.frailty_col,
         border = "grey30", lwd = 0.4, main = main, ...)
    add_frailty_legend(abs_max, col_seq, n_cols)
    return(invisible(data.frame(group = labs, mean = means, col = fill_cols)))
  }

  # ---- maps back-end -------------------------------------------------------
  if (!is.null(regions)) {
    if (!requireNamespace("maps", quietly = TRUE))
      stop("Package 'maps' is required for the maps back-end.", call. = FALSE)
    if (length(regions) != length(means))
      stop("regions must have the same length as group_levels.", call. = FALSE)

    world <- maps::map("world", plot = FALSE)
    maps::map("world", col = "grey90", fill = TRUE, border = "grey70",
              lwd = 0.3, main = main, ...)
    graphics::title(main = main)
    for (i in seq_along(regions)) {
      maps::map("world", regions = regions[i], col = fill_cols[i],
                fill = TRUE, border = "grey30", lwd = 0.5, add = TRUE)
    }
    add_frailty_legend(abs_max, col_seq, n_cols)
    return(invisible(data.frame(group = labs, mean = means, col = fill_cols)))
  }

  stop("Provide either 'regions' (maps back-end) or 'sf_map' + 'sf_key' (sf back-end).",
       call. = FALSE)
}

# -- Internal helpers --------------------------------------------------------

add_frailty_legend <- function(abs_max, col_seq, n_cols) {
  n_ticks <- 5
  tick_vals <- seq(-abs_max, abs_max, length.out = n_ticks)
  tick_labs <- sprintf("%.2f", tick_vals)

  usr <- graphics::par("usr")
  x0 <- usr[1] + (usr[2] - usr[1]) * 0.02
  x1 <- usr[1] + (usr[2] - usr[1]) * 0.04
  y0 <- usr[3] + (usr[4] - usr[3]) * 0.15
  y1 <- usr[3] + (usr[4] - usr[3]) * 0.85
  ys <- seq(y0, y1, length.out = n_cols + 1)

  for (k in seq_len(n_cols))
    graphics::rect(x0, ys[k], x1, ys[k + 1], col = col_seq[k], border = NA)
  graphics::rect(x0, y0, x1, y1, border = "grey30", lwd = 0.5)

  tick_y <- seq(y0, y1, length.out = n_ticks)
  graphics::segments(x1, tick_y, x1 + (x1 - x0) * 0.3, tick_y, col = "grey30")
  graphics::text(x1 + (x1 - x0) * 0.5, tick_y, tick_labs,
                 adj = 0, cex = 0.6, col = "grey30")
}

coef_labels <- function(terms_obj) {
  vars <- attr(terms_obj, "term.labels")
  c("(Intercept)", vars)
}
