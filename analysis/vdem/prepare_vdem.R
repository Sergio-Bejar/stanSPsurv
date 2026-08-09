#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(root, "analysis", "vdem", "output")
raw_dir <- file.path(root, "analysis", "vdem", "data-raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

needed <- c("countrycode", "maps", "survival")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

find_ert <- function(source = Sys.getenv("ERT_SOURCE", "")) {
  candidates <- c(
    source,
    file.path(raw_dir, "ERT"),
    file.path(raw_dir, "ERT_16.0")
  )
  candidates <- candidates[nzchar(candidates)]
  good <- candidates[file.exists(file.path(candidates, "data", "vdem.RData"))]
  if (length(good)) return(normalizePath(good[[1]], mustWork = TRUE))

  archive <- file.path(raw_dir, "ERT_16.0.tar.gz")
  url <- "https://github.com/vdeminstitute/ERT/raw/refs/heads/master/ERT_16.0.tar.gz"
  message("Downloading official ERT v16 source archive...")
  utils::download.file(url, archive, mode = "wb", quiet = FALSE)
  utils::untar(archive, exdir = raw_dir)
  unpacked <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
  good <- unpacked[file.exists(file.path(unpacked, "data", "vdem.RData"))]
  if (!length(good)) stop("Could not locate data/vdem.RData after unpacking ERT.")
  normalizePath(good[[1]], mustWork = TRUE)
}

zscore <- function(x) as.numeric(scale(x))

ert_root <- find_ert()
vdem_env <- new.env(parent = emptyenv())
load(file.path(ert_root, "data", "vdem.RData"), envir = vdem_env)
if (!exists("vdem", envir = vdem_env, inherits = FALSE)) stop("ERT vdem object not found.")
vdem <- get("vdem", envir = vdem_env)
ert <- utils::read.csv(file.path(ert_root, "inst", "ert.csv"), check.names = FALSE)

vars <- c(
  "country_name", "country_text_id", "country_id", "year", "COWcode",
  "v2x_regime", "v2x_jucon", "v2xlg_legcon", "v2x_freexp_altinf",
  "e_gdppc", "e_civil_war", "e_regionpol"
)
absent <- setdiff(vars, names(vdem))
if (length(absent)) stop("Required V-Dem variables missing: ", paste(absent, collapse = ", "))
d <- vdem[vdem$year >= 1900, vars]
d <- d[order(d$country_id, d$year), ]
d$democracy <- d$v2x_regime %in% c(2, 3)

# Lags are computed before restricting to democratic observations.
lag_within <- function(x, id) {
  out <- rep(NA_real_, length(x))
  same <- c(FALSE, id[-1] == id[-length(id)])
  out[same] <- x[which(same) - 1L]
  out
}
for (nm in c("v2x_jucon", "v2xlg_legcon", "v2x_freexp_altinf", "e_gdppc", "e_civil_war")) {
  d[[paste0("lag_", nm)]] <- lag_within(d[[nm]], d$country_id)
}

# A new spell starts at every transition from non-democracy (or missing) into democracy.
prev_dem <- ave(d$democracy, d$country_id, FUN = function(x) c(FALSE, head(x, -1)))
new_spell <- d$democracy & !prev_dem
d$spell_seq <- ave(as.integer(new_spell), d$country_id, FUN = cumsum)
dd <- d[d$democracy, ]
dd$spell_id <- paste(dd$country_text_id, dd$spell_seq, sep = "_")

make_spell <- function(x) {
  x <- x[order(x$year), ]
  id_rows <- which(d$country_id == x$country_id[1])
  last_global <- id_rows[d$year[id_rows] == max(x$year)][1]
  next_row <- if (last_global < nrow(d) && d$country_id[last_global + 1] == x$country_id[1]) last_global + 1 else NA_integer_
  breakdown <- !is.na(next_row) && !is.na(d$v2x_regime[next_row]) && d$v2x_regime[next_row] %in% c(0, 1)
  start <- x[1, ]
  data.frame(
    spell_id = start$spell_id,
    country_name = start$country_name,
    country_text_id = start$country_text_id,
    country_id = start$country_id,
    COWcode = start$COWcode,
    start_year = min(x$year),
    end_year = max(x$year),
    duration = nrow(x),
    breakdown = as.integer(breakdown),
    region = start$e_regionpol,
    judicial_constraints = start$lag_v2x_jucon,
    legislative_constraints = start$lag_v2xlg_legcon,
    free_expression = start$lag_v2x_freexp_altinf,
    log_gdppc = log(start$lag_e_gdppc),
    civil_war = start$lag_e_civil_war,
    stringsAsFactors = FALSE
  )
}
spells <- do.call(rbind, lapply(split(dd, dd$spell_id), make_spell))
row.names(spells) <- NULL

# Democratic years accumulated before the current spell.
spells <- spells[order(spells$country_id, spells$start_year), ]
spells$prior_dem_years <- ave(spells$duration, spells$country_id,
                             FUN = function(x) c(0, head(cumsum(x), -1)))

# Attach ERT's classification at the spell onset for validation and description.
ert_keep <- ert[, intersect(c("country_id", "year", "dem_ep", "dem_ep_id",
                              "dem_ep_outcome", "dem_ep_outcome_agg"), names(ert))]
names(ert_keep)[names(ert_keep) == "year"] <- "start_year"
spells <- merge(spells, ert_keep, by = c("country_id", "start_year"), all.x = TRUE)

model_vars <- c("judicial_constraints", "legislative_constraints", "free_expression",
                "log_gdppc", "civil_war", "prior_dem_years")
complete <- stats::complete.cases(spells[, model_vars]) & spells$duration > 0
analysis <- spells[complete, ]
for (nm in model_vars) analysis[[paste0(nm, "_z")]] <- zscore(analysis[[nm]])

# Map country labels to polygon centroids and create a connected geographic graph.
world <- maps::map("world", plot = FALSE, fill = TRUE)
poly_name <- sub(":.*$", "", world$names)
# maps stores polygons separated by NA; compute each named polygon directly for robustness.
cent <- do.call(rbind, lapply(unique(poly_name), function(nm) {
  p <- maps::map("world", regions = nm, plot = FALSE, fill = TRUE)
  data.frame(map_name = nm, lon = mean(p$x, na.rm = TRUE), lat = mean(p$y, na.rm = TRUE))
}))
cent$iso3c <- countrycode::countrycode(cent$map_name, "country.name", "iso3c", warn = FALSE)
cent <- aggregate(cbind(lon, lat) ~ iso3c, cent, mean)
analysis$iso3c <- countrycode::countrycode(analysis$country_name, "country.name", "iso3c", warn = FALSE)
areas <- sort(unique(analysis$country_text_id[!is.na(analysis$iso3c)]))
loc <- unique(analysis[match(areas, analysis$country_text_id), c("country_text_id", "iso3c")])
loc <- merge(loc, cent, by = "iso3c", all.x = TRUE)
loc <- loc[match(areas, loc$country_text_id), ]
has_location <- stats::complete.cases(loc[, c("lon", "lat")])
areas <- areas[has_location]
loc <- loc[has_location, ]
analysis <- analysis[analysis$country_text_id %in% areas, ]

great_circle <- function(lon, lat) {
  rad <- pi / 180
  lon <- lon * rad; lat <- lat * rad
  cosine <- outer(sin(lat), sin(lat)) + outer(cos(lat), cos(lat)) *
    cos(outer(lon, lon, "-"))
  cosine[cosine > 1] <- 1
  cosine[cosine < -1] <- -1
  acos(cosine)
}
dist <- great_circle(loc$lon, loc$lat)
diag(dist) <- Inf
k <- min(4L, nrow(loc) - 1L)
A <- matrix(0L, nrow(loc), nrow(loc), dimnames = list(areas, areas))
for (i in seq_len(nrow(loc))) A[i, order(dist[i, ])[seq_len(k)]] <- 1L
A <- 1L * ((A + t(A)) > 0)
diag(A) <- 0L

# Ordering matters: factor levels in stanSPsurv are alphabetical.
stopifnot(identical(rownames(A), sort(unique(analysis$country_text_id))))

metadata <- list(
  data_release = "V-Dem 16 / ERT 16.0",
  years = range(analysis$start_year),
  outcome = "transition from v2x_regime 2/3 to 0/1",
  spatial_graph = paste0(k, "-nearest polygon-centroid neighbours, symmetrized"),
  created = as.character(Sys.Date())
)
saveRDS(list(data = analysis, adjacency = A, locations = loc, metadata = metadata),
        file.path(out_dir, "vdem_application.rds"))

summary <- data.frame(
  statistic = c("spells_all", "spells_analysis", "breakdowns_analysis", "countries_analysis",
                "first_start_year", "last_start_year", "adjacency_edges", "isolated_areas"),
  value = c(nrow(spells), nrow(analysis), sum(analysis$breakdown),
            length(unique(analysis$country_text_id)), min(analysis$start_year),
            max(analysis$start_year), sum(A) / 2, sum(rowSums(A) == 0))
)
utils::write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
message("Prepared ", nrow(analysis), " spells (", sum(analysis$breakdown),
        " breakdowns) across ", length(unique(analysis$country_text_id)), " countries.")
