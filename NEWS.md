/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
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
  spells, including a geographic-proximity graph and ready-to-run pooled,
  exchangeable, and spatial specifications.
- Recorded the first full pooled Weibull application run and its sampler
  diagnostics.

- Created the initial R package structure.
- Added candidate pooled, exchangeable, and spatial split-population survival models.
- Added Weibull and log-logistic baselines through a shared Stan backend.
- Added coefficient, prediction, and frailty helpers.
- Added initial input-validation tests and package branding.
- Compiled the shared backend with CmdStan 2.39.0 and ran all six candidate
  model combinations on simulated data.
- Fixed an improper auxiliary parameter block in the pooled specification.
- Added a reproducible validation script and recorded convergence diagnostics.
- Achieved `R CMD check --no-manual` status OK on R 4.2.1 for macOS.
- Added counting-process `Surv(t0, time, event)` support for time-varying
  country-year covariates and repeated episode rows.
- Added distinct split-stage and survival-stage frailties with separate CAR
  processes in spatial models.
- Extended `frailty()` with a `stage` argument and added start-stop tests.
