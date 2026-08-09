# stanSPsurv

<p align="center">
  <img src="man/figures/logo.png" alt="stanSPsurv logo" width="260">
</p>

`stanSPsurv` is a development-stage R package for Bayesian
split-population survival models fitted with Stan. It provides a common API
for pooled, exchangeable-frailty, and spatial CAR specifications with
Weibull or log-logistic event-time distributions.

> **Validation status:** all six candidate specifications compile and complete
> short simulation runs without divergences, and `R CMD check --no-manual`
> passes with `Status: OK`. The dual-CAR spatial smoke tests reached maximum
> tree depth in 32/600 (Weibull) and 119/600 (log-logistic) transitions. Larger simulation
> studies and comparison with `BayesSPsurv` remain required before substantive
> use.

> **Walter benchmark:** the real-data pooled Weibull run produced 394/2,000
> divergences and the exchangeable run produced 542/2,000. The official 800 km
> graph contains 16 isolated countries, and the spatial smoke run stalled after
> entering sampling. These results identify a package-level parameterization
> problem; do not interpret real-data spatial estimates yet. See
> [`validation/walter`](validation/walter).

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

Both `Surv(time, event)` and counting-process
`Surv(t0, time, event)` responses are supported. Exchangeable and spatial
models estimate separate frailties for the split and survival stages; use
`frailty(fit, stage = "split")` or `frailty(fit, stage = "survival")`.

## Reproducible validation

Run the six-model smoke/recovery study from the repository root:

```sh
Rscript validation/run_models.R
```

The concise results are stored in `validation/results.csv`; raw CmdStan output
is intentionally excluded from version control.

## Political-science application

A reproducible V-Dem 16 / ERT 16 pipeline for democratic-breakdown spells is
in [`analysis/vdem`](analysis/vdem). It prepares pooled, country-frailty, and
spatial-proximity specifications while keeping raw source data and fitted
objects out of version control.

## Authors

Sergio Béjar (first author and maintainer) and Bumba Mukherjee (coauthor).
