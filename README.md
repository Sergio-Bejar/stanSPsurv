# stanSPsurv

<p align="center">
  <img src="man/figures/logo.png" alt="stanSPsurv logo" width="260">
</p>

`stanSPsurv` is a development-stage R package for Bayesian
split-population survival models fitted with Stan. It provides a common API
for pooled, exchangeable-frailty, and spatial CAR specifications with
Weibull or log-logistic event-time distributions.

> **Validation status:** the package architecture and candidate likelihood are
> implemented, but numerical recovery, comparison with `BayesSPsurv`, and a
> complete `R CMD check` remain required before substantive use.

## Installation

The package is currently private and under active development. It requires R
4.1 or newer, `cmdstanr`, and a working CmdStan installation.

```r
cmdstanr::install_cmdstan()
```

## Example

```r
library(survival)
library(stanSPsurv)

fit <- stan_pooled_spsurv(
  Surv(time, status) ~ age + sex,
  split_formula = ~ age + sex,
  data = lung,
  form = "weibull",
  chains = 4,
  parallel_chains = 4
)

summary(fit)
coef(fit, component = "survival")
```

## Model families

- `stan_pooled_spsurv()`
- `stan_exchange_spsurv()`
- `stan_spatial_spsurv()`

Each supports `form = "weibull"` or `form = "loglogistic"`.

## Authors

Sergio Béjar (first author and maintainer) and Bumba Mukherjee (coauthor).

