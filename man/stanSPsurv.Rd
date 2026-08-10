/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
\name{stanSPsurv-models}
\alias{stan_pooled_spsurv}
\alias{stan_exchange_spsurv}
\alias{stan_spatial_spsurv}
\title{Fit Bayesian split-population survival models with Stan}
\description{
Fit pooled, exchangeable-frailty, or spatial CAR split-population survival
models with Weibull or log-logistic event-time distributions.
}
\usage{
stan_pooled_spsurv(formula, split_formula, data, spell = NULL,
  form = c("weibull", "loglogistic"), ...)

stan_exchange_spsurv(formula, split_formula, data, group, spell = NULL,
  form = c("weibull", "loglogistic"), ...)

stan_spatial_spsurv(formula, split_formula, data, area, adjacency, spell = NULL,
  form = c("weibull", "loglogistic"), ...)
}
\arguments{
  \item{formula}{Survival-stage formula with a right-censored
  \code{survival::Surv(time, event)} or
  \code{survival::Surv(t0, time, event)} response.}
  \item{split_formula}{Formula for the latent susceptibility probability.}
  \item{data}{A data frame.}
  \item{spell}{Episode identifier. Required with counting-process responses.
  Split-stage covariates must be constant within an episode.}
  \item{form}{The event-time distribution.}
  \item{group}{Group identifier as a column name or vector.}
  \item{area}{Spatial unit identifier as a column name or vector.}
  \item{adjacency}{Symmetric binary adjacency matrix. Disconnected components
  and isolated areas are supported by the proper-CAR mixture.}
  \item{...}{Sampling arguments passed to CmdStan.}
}
\value{An object of class \code{stanSPsurv}.}
