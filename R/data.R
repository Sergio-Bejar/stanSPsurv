prepare_spsurv_data <- function(formula, split_formula, data, group = NULL,
                                adjacency = NULL, spell = NULL) {
  mf_y <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf_y)
  if (!inherits(y, "Surv") || !ncol(y) %in% c(2L, 3L))
    stop("formula must use Surv(time, event) or Surv(t0, time, event).", call. = FALSE)
  mf_z <- stats::model.frame(split_formula, data = data, na.action = stats::na.pass)
  keep <- if (ncol(mf_z) > 0L) stats::complete.cases(mf_y, mf_z) else stats::complete.cases(mf_y)
  if (!is.null(group)) keep <- keep & !is.na(group)
  if (!is.null(spell)) keep <- keep & !is.na(spell)
  if (!all(keep)) {
    mf_y <- mf_y[keep, , drop = FALSE]
    mf_z <- mf_z[keep, , drop = FALSE]
    y <- y[keep, , drop = FALSE]
    if (!is.null(group)) group <- group[keep]
    if (!is.null(spell)) spell <- spell[keep]
  }
  X <- stats::model.matrix(stats::terms(mf_y), mf_y)
  Z <- stats::model.matrix(stats::terms(mf_z), mf_z)
  if (ncol(y) == 2L) {
    t0 <- rep(0, nrow(y)); time <- y[, 1]; event <- y[, 2]
  } else {
    if (is.null(spell))
      stop("spell is required with Surv(t0, time, event).", call. = FALSE)
    t0 <- y[, 1]; time <- y[, 2]; event <- y[, 3]
  }
  if (any(t0 < 0) || any(time <= t0))
    stop("Intervals must satisfy 0 <= t0 < time.", call. = FALSE)
  if (qr(X)$rank < ncol(X) || qr(Z)$rank < ncol(Z))
    stop("Design matrices must be full rank.", call. = FALSE)
  if (is.null(group)) {
    group_id <- rep.int(1L, nrow(X)); group_levels <- NULL; J <- 1L
  } else {
    f <- factor(group); group_id <- as.integer(f); group_levels <- levels(f); J <- nlevels(f)
  }
  if (is.null(spell)) spell <- seq_len(nrow(X))
  spell_factor <- factor(spell)
  order_rows <- order(as.integer(spell_factor), t0, time)
  X <- X[order_rows, , drop = FALSE]
  Z <- Z[order_rows, , drop = FALSE]
  t0 <- t0[order_rows]; time <- time[order_rows]; event <- event[order_rows]
  group_id <- group_id[order_rows]
  spell_factor <- droplevels(spell_factor[order_rows])
  spell_id <- as.integer(spell_factor)
  rows <- split(seq_along(spell_id), spell_id)
  bad_event <- vapply(rows, function(ii) sum(event[ii]) > 1L ||
                        (any(event[ii] == 1L) && utils::tail(event[ii], 1) != 1L), logical(1))
  if (any(bad_event)) stop("Each spell may have one event, only in its final interval.", call. = FALSE)
  bad_group <- vapply(rows, function(ii) length(unique(group_id[ii])) != 1L, logical(1))
  if (any(bad_group)) stop("group/area must be constant within each spell.", call. = FALSE)
  bad_z <- vapply(rows, function(ii) any(vapply(seq_len(ncol(Z)), function(j)
    length(unique(Z[ii, j])) != 1L, logical(1))), logical(1))
  if (any(bad_z)) stop("split_formula covariates must be constant within each spell.", call. = FALSE)
  first <- vapply(rows, `[`, integer(1), 1L)
  last <- vapply(rows, function(ii) utils::tail(ii, 1), integer(1))
  Z_spell <- Z[first, , drop = FALSE]
  spell_group_id <- group_id[first]
  spell_event <- as.integer(vapply(rows, function(ii) any(event[ii] == 1L), logical(1)))
  edges <- adjacency_edges(adjacency)
  graph <- adjacency_graph(adjacency)
  list(stan_data = list(N = nrow(X), M = length(rows), K = ncol(X), Q = ncol(Z),
       t0 = as.vector(t0), time = as.vector(time), event = as.integer(event), X = X, Z = Z,
       J = J, group_id = group_id, N_edges = length(edges$node1),
       node1 = edges$node1, node2 = edges$node2,
       spell_start = first, spell_end = last, spell_event = spell_event,
       Z_spell = Z_spell, spell_group_id = spell_group_id,
       component_id = graph$component_id, N_components = graph$N_components,
       component_size = graph$component_size, icar_scale = graph$scale),
       terms = stats::terms(mf_y), split_terms = stats::terms(mf_z),
       group_levels = group_levels)
}

validate_adjacency <- function(A, J) {
  if (!is.matrix(A) || any(dim(A) != c(J, J)))
    stop("adjacency must be a J by J matrix.", call. = FALSE)
  if (!isTRUE(all.equal(A, t(A))) || any(diag(A) != 0) || any(!A %in% c(0, 1)))
    stop("adjacency must be symmetric, binary, and have a zero diagonal.", call. = FALSE)
  invisible(TRUE)
}

adjacency_edges <- function(A) {
  if (is.null(A)) return(list(node1 = integer(), node2 = integer()))
  idx <- which(upper.tri(A) & A == 1, arr.ind = TRUE)
  list(node1 = as.integer(idx[, 1]), node2 = as.integer(idx[, 2]))
}

adjacency_graph <- function(A) {
  if (is.null(A)) return(list(component_id = 1L, N_components = 1L,
                              component_size = 1L, scale = 1))
  J <- nrow(A); component <- integer(J); current <- 0L
  for (root in seq_len(J)) if (component[root] == 0L) {
    current <- current + 1L; queue <- root; component[root] <- current
    while (length(queue)) {
      i <- queue[1]; queue <- queue[-1]
      nbr <- which(A[i, ] == 1L & component == 0L)
      if (length(nbr)) { component[nbr] <- current; queue <- c(queue, nbr) }
    }
  }
  sizes <- tabulate(component, current)
  scales <- rep(1, J)
  for (cc in seq_len(current)) {
    idx <- which(component == cc)
    if (length(idx) > 1L) {
      Q <- diag(rowSums(A[idx, idx, drop = FALSE])) - A[idx, idx, drop = FALSE]
      ev <- eigen(Q, symmetric = TRUE)
      keep <- ev$values > sqrt(.Machine$double.eps)
      Qplus <- ev$vectors[, keep, drop = FALSE] %*%
        diag(1 / ev$values[keep], nrow = sum(keep)) %*%
        t(ev$vectors[, keep, drop = FALSE])
      scales[idx] <- exp(mean(log(diag(Qplus))))
    }
  }
  list(component_id = component, N_components = current,
       component_size = sizes, scale = scales)
}
