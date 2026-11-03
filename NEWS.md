# stanSPsurv 0.1.0

## CRAN submission

- First public release. `R CMD check --as-cran` passes with Status: OK.
- License changed from CC0 to MIT.

## Validation

- Added 100-replication parameter recovery study for all six specifications
  (3 structures × 2 baselines, N = 500). Zero divergences across 600 runs,
  Rhat < 1.001, 81–97 % empirical coverage of 90 % credible intervals.
- Added 50-replication head-to-head comparison with `BayesSPsurv`:
  comparable bias, 40–1 000× higher ESS per second in stanSPsurv.
- Documented a case-sensitivity bug in BayesSPsurv's MCMC engine that
  prevents the Weibull shape parameter from being sampled.

## Plotting

- Added `plot.stanSPsurv()` for coefficient forest plots.
- Added `plot_survival()` for survival curves with credible bands.
- Added `plot_susceptibility()` for susceptibility probability curves.
- Added `plot_frailty()` for frailty caterpillar plots.
- Added `plot_choropleth()` for spatial frailty maps with `maps` and
  `sf` back-ends.

## Data and documentation

- Added `conflict_spells` simulated dataset for examples and testing.
- Added introductory vignette covering all model families and plotting.
- Full roxygen documentation for all 13 exported functions and methods.
- Added CITATION file.

## Bug fixes

- Fixed `complete.cases()` failure when split formula is intercept-only.

# stanSPsurv 0.0.3.9000

- Split the pooled, exchangeable, and spatial implementations into separate
  Stan programs so inactive frailty parameters are never sampled.
- Corrected counting-process likelihoods to mix susceptible and immune
  populations once per episode through the new `spell` argument.
- Reparameterized positive scales on the log scale and added regularizing
  priors for weakly identified split-population models.
- Replaced the unscaled dual ICAR with a scaled proper-CAR mixture that handles
  disconnected components and isolated areas.
- Re-ran the Walter (2015) benchmark: pooled, exchangeable, and spatial models
  completed with zero divergences and zero maximum-treedepth transitions.

# stanSPsurv 0.0.2.9000

- Added a reproducible V-Dem 16 / ERT 16 application for democratic-breakdown
  spells.
- Created the initial R package structure with pooled, exchangeable, and
  spatial split-population survival models.
- Added Weibull and log-logistic baselines, coefficient/prediction/frailty
  helpers, counting-process support, and input-validation tests.
