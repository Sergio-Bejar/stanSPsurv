#' @keywords internal
"_PACKAGE"

#' Simulated post-conflict survival data
#'
#' A simulated dataset of 300 post-conflict episodes across 20 countries,
#' designed to illustrate split-population survival models. Countries are
#' grouped into 4 geographic regions with a chain adjacency structure.
#'
#' @format A data frame with 300 rows and 8 columns:
#' \describe{
#'   \item{time}{Duration until war recurrence or censoring.}
#'   \item{event}{Binary indicator: 1 if war recurred, 0 if censored.}
#'   \item{victory}{Binary: 1 if the previous war ended in decisive victory.}
#'   \item{gdp}{Log GDP per capita (standardized).}
#'   \item{peacekeeping}{Binary: 1 if UN peacekeepers were deployed.}
#'   \item{country}{Factor with 20 country labels.}
#'   \item{region}{Factor with 4 region labels.}
#'   \item{country_id}{Integer country identifier (1--20).}
#' }
#'
#' @details
#' The data are entirely simulated with known parameter values
#' (beta = c(-2.0, -0.5, 0.3, -0.8), gamma = c(0.8, -0.6),
#' shape = 1.2) and are intended for package examples and testing.
#' They do not represent real conflict data.
#'
#' @examples
#' data(conflict_spells)
#' head(conflict_spells)
#'
#' @source Simulated for package illustration.
"conflict_spells"
