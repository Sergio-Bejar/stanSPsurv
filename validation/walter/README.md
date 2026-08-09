/Users/sergiobejar/.rvm/scripts/rvm:29: operation not permitted: ps
# Walter (2015) benchmark

This benchmark reproduces the data preparation and model formula documented by
`BayesSPsurv` for Barbara Walter's post-war peace data. The input has 1,237
country-years, 61 conflict spells, 46 countries, and 24 recurrent-war events.
The published 800 km capital-distance graph has 27 undirected edges and 16
isolated countries.

The Stan runs use four chains, `adapt_delta = 0.97`, and
`max_treedepth = 13`. Pooled and exchangeable runs contain 500 warm-up and 500
sampling iterations per chain. The spatial diagnostic was shortened to 250 +
250 and stopped after all chains entered sampling but failed to complete the
next 50-iteration reporting block within several minutes.

## Results

| Specification | Divergences | Max treedepth | R-hat | Assessment |
|---|---:|---:|---:|---|
| Pooled Weibull | 394/2,000 | 0 | <= 1.01 | Failed: divergent geometry without CAR |
| Exchangeable Weibull | 542/2,000 | 0 | <= 1.01 | Failed: frailties worsen divergences |
| Spatial, official 800 km | not completed | stalled | not available | Failed computational smoke test |

The pooled failure establishes that Walter-specific spatial adjacency is not
the root cause. In the current Stan program, `sigma_v`, `sigma_w`, `v_raw`, and
`w_raw` are sampled even for the pooled model although they do not enter its
likelihood. The constrained scale parameters can create boundary geometry with
no compensating information. For spatial models, the dual unscaled ICAR and
the 16 isolated nodes add further identification problems.

Before interpreting Walter coefficients, the package should use
structure-specific Stan programs (or conditionally sized parameter blocks), a
non-centered exchangeable parameterization with regularizing priors, and a
scaled BYM2 component-aware spatial prior. The benchmark should then be rerun
against the original `BayesSPsurv` estimates.

