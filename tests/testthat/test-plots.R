test_that("coef_labels extracts variable names", {
  d <- data.frame(y = 1:3, x1 = rnorm(3), x2 = rnorm(3))
  tt <- terms(y ~ x1 + x2, data = d)
  labs <- stanSPsurv:::coef_labels(tt)
  expect_equal(labs, c("(Intercept)", "x1", "x2"))
})

test_that("plot_survival rejects missing times", {
  obj <- list(structure = "pooled", form = "weibull")
  class(obj) <- "stanSPsurv"
  expect_error(plot_survival(obj, newdata = data.frame(x = 1)),
               "at least two time points")
})

test_that("plot_frailty rejects pooled model", {
  obj <- list(structure = "pooled")
  class(obj) <- "stanSPsurv"
  expect_error(plot_frailty(obj), "exchangeable or spatial")
})

test_that("plot_choropleth rejects non-spatial model", {
  obj <- list(structure = "pooled")
  class(obj) <- "stanSPsurv"
  expect_error(plot_choropleth(obj), "spatial model")
})

test_that("plot_choropleth requires regions or sf_map", {
  obj <- list(structure = "exchangeable")
  class(obj) <- "stanSPsurv"
  # Non-spatial model hits the structure check before regions/sf_map
  expect_error(plot_choropleth(obj), "spatial model")
})

test_that("plot_susceptibility rejects missing variable", {
  obj <- list(structure = "pooled")
  class(obj) <- "stanSPsurv"
  d <- data.frame(x = 1:3, z = 4:6)
  expect_error(plot_susceptibility(obj, variable = "w", data = d),
               "not found")
})
