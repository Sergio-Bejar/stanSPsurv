/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
structure <- if (length(args)) args[[1]] else "pooled"
graph <- if (length(args) > 1) args[[2]] else "official800"
stopifnot(structure %in% c("pooled", "exchangeable", "spatial"))
stopifnot(graph %in% c("official800", "connected"))
reference <- Sys.getenv("BAYESSPSURV_REFERENCE", "/private/tmp/bayesspsurv-reference")
e <- new.env(parent = emptyenv())
load(file.path(reference, "data", "Walter_2015_JCR.rda"), envir = e)
load(file.path(reference, "data", "capdist.rda"), envir = e)
d <- spduration::add_duration(e$Walter_2015_JCR, "renewed_war",
                              unitID = "id", tID = "year",
                              freq = "year", ongoing = FALSE)
cc <- sort(unique(d$ccode))
A <- matrix(0L, length(cc), length(cc), dimnames = list(cc, cc))
dist <- e$capdist[e$capdist$numa %in% cc & e$capdist$numb %in% cc,
                  c("numa", "numb", "kmdist")]
near <- dist[dist$kmdist <= 800, ]
A[cbind(match(near$numa, cc), match(near$numb, cc))] <- 1L
A <- 1L * ((A + t(A)) > 0)
diag(A) <- 0L
if (graph == "connected") {
  while (any(rowSums(A) == 0)) {
    i <- which(rowSums(A) == 0)[1]
    choices <- dist[dist$numa == cc[i] & dist$numb %in% cc, ]
    choices <- choices[order(choices$kmdist), ]
    j <- match(choices$numb[1], cc)
    A[i, j] <- A[j, i] <- 1L
  }
}
f_duration <- survival::Surv(t.0, duration, failure) ~
  victory + comprehensive + lgdpl + unpko
f_split <- ~ lgdpl
warmup <- as.integer(Sys.getenv("WALTER_WARMUP", "500"))
sampling <- as.integer(Sys.getenv("WALTER_SAMPLING", "500"))
common <- list(form = "weibull", chains = 4, parallel_chains = 4,
               iter_warmup = warmup, iter_sampling = sampling, seed = 2015,
               adapt_delta = 0.97, max_treedepth = 13, refresh = 100)
if (structure == "pooled") {
  fit <- do.call(stanSPsurv::stan_pooled_spsurv,
                 c(list(f_duration, f_split, d), common))
} else if (structure == "exchangeable") {
  fit <- do.call(stanSPsurv::stan_exchange_spsurv,
                 c(list(f_duration, f_split, d, group = d$ccode), common))
} else {
  fit <- do.call(stanSPsurv::stan_spatial_spsurv,
                 c(list(f_duration, f_split, d, area = d$ccode,
                        adjacency = A), common))
}
out <- file.path("validation", "walter", "output")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stem <- paste(structure, graph, "weibull", sep = "_")
s <- summary(fit)
diagnostics <- fit$fit$diagnostic_summary()
diag_table <- data.frame(chain = seq_along(diagnostics$num_divergent),
                         divergent = diagnostics$num_divergent,
                         max_treedepth = diagnostics$num_max_treedepth,
                         ebfmi = diagnostics$ebfmi)
utils::write.csv(s, file.path(out, paste0(stem, "_summary.csv")), row.names = FALSE)
utils::write.csv(diag_table, file.path(out, paste0(stem, "_diagnostics.csv")), row.names = FALSE)
utils::write.csv(data.frame(rows = nrow(d), spells = length(unique(d$id)),
                            countries = length(cc), events = sum(d$failure),
                            edges = sum(A) / 2, isolates = sum(rowSums(A) == 0)),
                 file.path(out, paste0(stem, "_data.csv")), row.names = FALSE)
print(s)
print(diag_table)
