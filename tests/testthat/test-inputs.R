test_that("adjacency validation rejects malformed matrices", {
  expect_error(stanSPsurv:::validate_adjacency(diag(2), 2), "zero diagonal")
  expect_silent(stanSPsurv:::validate_adjacency(matrix(c(0,1,1,0), 2), 2))
})

test_that("start-stop survival intervals are prepared", {
  skip_if_not_installed("survival")
  d <- data.frame(t0 = c(0, 1), t = c(1, 2), event = c(0, 1), x = c(-1, 1))
  out <- stanSPsurv:::prepare_spsurv_data(
    survival::Surv(t0, t, event) ~ x, ~ 1, d, spell = c(1, 1)
  )
  expect_equal(out$stan_data$t0, d$t0)
  expect_equal(out$stan_data$time, d$t)
})

test_that("start-stop responses require spell identifiers", {
  skip_if_not_installed("survival")
  d <- data.frame(t0 = c(0, 1), t = c(1, 2), event = c(0, 1), x = c(-1, 1))
  expect_error(stanSPsurv:::prepare_spsurv_data(
    survival::Surv(t0, t, event) ~ x, ~ 1, d), "spell is required")
})

test_that("split covariates are spell-constant", {
  skip_if_not_installed("survival")
  d <- data.frame(t0 = c(0, 1), t = c(1, 2), event = c(0, 1), x = c(-1, 1))
  expect_error(stanSPsurv:::prepare_spsurv_data(
    survival::Surv(t0, t, event) ~ x, ~ x, d, spell = c(1, 1)),
    "constant within each spell")
})

test_that("intercept-only split formula works", {
  skip_if_not_installed("survival")
  d <- data.frame(time = c(1, 2, 3), event = c(0, 1, 0), x = c(-1, 0, 1))
  out <- stanSPsurv:::prepare_spsurv_data(
    survival::Surv(time, event) ~ x, ~ 1, d
  )
  expect_equal(out$stan_data$N, 3L)
  expect_equal(ncol(out$stan_data$Z_spell), 1L)
})
