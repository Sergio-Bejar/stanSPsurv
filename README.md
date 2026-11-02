# stanSPsurv

<p align="center">
  <img src="man/figures/logo.png" alt="stanSPsurv logo" width="260">
</p>

`stanSPsurv` is a development-stage R package for Bayesian
split-population survival models fitted with Stan. It provides a common API
for pooled, exchangeable-frailty, and scaled proper-CAR specifications with
Weibull or log-logistic event-time distributions.

> **Validation status:** all six specifications pass a 100-replication
> parameter recovery study (N = 500) with zero divergences, Rhat < 1.001,
> and 81–97 % empirical coverage of 90 % credible intervals. A formal
> 50-replication comparison with `BayesSPsurv` shows comparable bias when
> the Weibull shape is within BayesSPsurv's support (< 1), with stanSPsurv
> delivering 40–1 000× higher ESS per second. `R CMD check --no-manual`
> passes with `Status: OK`. See [Validation](#validation) below.

> **Walter benchmark:** pooled, exchangeable, and 800 km spatial runs
> completed with zero divergences and zero maximum-treedepth transitions
> after introducing episode-level mixtures, separate Stan programs,
> log-scale variance parameters, and a proper-CAR mixture. See
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
`Surv(t0, time, event)` responses are supported. Counting-process models require
the `spell` argument and evaluate latent susceptibility once per episode.
Exchangeable and spatial
models estimate separate frailties for the split and survival stages; use
`frailty(fit, stage = "split")` or `frailty(fit, stage = "survival")`.

## Plotting

Five plotting functions are included, all using base R graphics.
The panel below shows all four on the Walter (2015) post-war peace data
(spatial Weibull model, 46 countries, 800 km adjacency graph):

<p align="center">
  <img src="man/figures/walter_panel.png" alt="Walter panel" width="900">
</p>

### Coefficient plot

```r
plot(fit)
```

### Survival curves

```r
plot_survival(fit, newdata = data.frame(victory = c(0, 1), comprehensive = c(0, 0),
                                        lgdpl_c = 0, unpko = 0),
              times = seq(0.5, 30, length.out = 60),
              labels = c("Baseline", "Victory"))
```

### Susceptibility and frailty

```r
plot_susceptibility(fit, variable = "lgdpl_split_c", data = d)
plot_frailty(fit, stage = "survival")
```

### Choropleth map of spatial frailties

```r
plot_choropleth(fit, stage = "survival", regions = country_names)
```

<p align="center">
  <img src="man/figures/choropleth_example.png" alt="Choropleth" width="740">
</p>

The `maps` back-end requires the `maps` package (in `Suggests`).
An `sf` back-end is also available for arbitrary spatial geometries:

```r
plot_choropleth(fit, stage = "survival", sf_map = my_sf, sf_key = "iso3")
```

The `maps` back-end requires the `maps` package (in `Suggests`).
An `sf` back-end is also available for arbitrary spatial geometries:

```r
plot_choropleth(fit, stage = "survival", sf_map = my_sf, sf_key = "iso3")
```

## Validation

### Parameter recovery (simulation study)

The full study runs 100 replications of each of the six model specifications
(3 structures × 2 baselines) on data simulated with known parameters
(N = 500 per replication):

```sh
Rscript validation/simulation_study.R          # defaults: --reps 100 --N 500
Rscript validation/simulation_study.R --reps 5 # quick smoke test
```

Results from the latest run (`validation/simulation_output/summary.csv`):

| Structure    | Baseline     | Max |bias| | Max RMSE | Coverage 90% | Divergences | Median sec |
|--------------|--------------|-------------|----------|--------------|-------------|------------|
| Pooled       | Weibull      | 0.014       | 0.110    | 87–93 %      | 0           | 4.6        |
| Pooled       | Log-logistic | 0.043       | 0.163    | 81–93 %      | 0           | 2.8        |
| Exchangeable | Weibull      | 0.025       | 0.135    | 86–90 %      | 0           | 9.7        |
| Exchangeable | Log-logistic | 0.069       | 0.181    | 90–95 %      | 0           | 5.8        |
| Spatial      | Weibull      | 0.031       | 0.133    | 85–97 %      | 0           | 236        |
| Spatial      | Log-logistic | 0.049       | 0.173    | 86–96 %      | 0           | 242        |

### Comparison with BayesSPsurv

A head-to-head comparison fits the same simulated data with both packages
(50 replications, N = 400, pooled and exchangeable Weibull with shape = 0.75):

```sh
Rscript validation/compare_bayesspsurv.R          # defaults: --reps 50 --N 400
Rscript validation/compare_bayesspsurv.R --reps 5 # quick smoke test
```

| Metric                    | stanSPsurv  | BayesSPsurv |
|---------------------------|-------------|-------------|
| Bias (pooled)             | comparable  | comparable  |
| Bias (exchangeable gamma) | 0.09        | 0.21        |
| Coverage (exch. gamma)    | 96 %        | 84 %        |
| ESS / sec (pooled)        | 487–582     | 13          |
| ESS / sec (exchangeable)  | 304–468     | 0.3–3       |
| Wall time (pooled)        | 3 sec       | 152 sec     |
| Wall time (exchangeable)  | 6 sec       | 239 sec     |

Two limitations of BayesSPsurv affect the comparison. First, its Weibull shape
parameter `rho` is constrained to `rho >= 0.01` with a Gamma(1,1) prior, so
`shape = exp(-rho) <= 0.99` — only decreasing hazards can be represented.
Second, a case-sensitivity bug in the MCMC engine (`form %in% "Weibull"` vs
the lowercase `"weibull"` passed by `pooledSPsurv()`) prevents `rho` from
being sampled at all unless patched locally.

### Smoke test (six-model)

The original quick smoke test is also available:

```sh
Rscript validation/run_models.R
```

Results in `validation/results.csv`.

## Political-science application

A reproducible V-Dem 16 / ERT 16 pipeline for democratic-breakdown spells is
in [`analysis/vdem`](analysis/vdem). It prepares pooled, country-frailty, and
spatial-proximity specifications while keeping raw source data and fitted
objects out of version control.

## Authors

Sergio Béjar (first author and maintainer) and Bumba Mukherjee (coauthor).
