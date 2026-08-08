prepare_spsurv_data <- function(formula, split_formula, data, group = NULL,
                                adjacency = NULL) {
  mf_y <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf_y)
  if (!inherits(y, "Surv") || ncol(y) != 2L)
    stop("formula must have a right-censored Surv(time, event) response.", call. = FALSE)
  mf_z <- stats::model.frame(split_formula, data = data, na.action = stats::na.pass)
  keep <- stats::complete.cases(mf_y, mf_z)
  if (!is.null(group)) keep <- keep & !is.na(group)
  if (!all(keep)) {
    mf_y <- mf_y[keep, , drop = FALSE]
    mf_z <- mf_z[keep, , drop = FALSE]
    y <- y[keep, , drop = FALSE]
    if (!is.null(group)) group <- group[keep]
  }
  X <- stats::model.matrix(stats::terms(mf_y), mf_y)
  Z <- stats::model.matrix(stats::terms(mf_z), mf_z)
  if (any(y[, 1] <= 0)) stop("All event/censoring times must be positive.", call. = FALSE)
  if (qr(X)$rank < ncol(X) || qr(Z)$rank < ncol(Z))
    stop("Design matrices must be full rank.", call. = FALSE)
  if (is.null(group)) {
    group_id <- rep.int(1L, nrow(X)); group_levels <- NULL; J <- 1L
  } else {
    f <- factor(group); group_id <- as.integer(f); group_levels <- levels(f); J <- nlevels(f)
  }
  edges <- adjacency_edges(adjacency)
  list(stan_data = list(N = nrow(X), K = ncol(X), Q = ncol(Z),
       time = as.vector(y[, 1]), event = as.integer(y[, 2]), X = X, Z = Z,
       J = J, group_id = group_id, N_edges = length(edges$node1),
       node1 = edges$node1, node2 = edges$node2),
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

