test_that("two-column Surv is prepared correctly", {
  skip_if_not_installed("survival")
  d <- data.frame(time = c(1, 5, 10), event = c(1, 0, 1),
                  x = c(-1, 0, 1), z = c(0.5, -0.5, 0))
  out <- stanSPsurv:::prepare_spsurv_data(
    survival::Surv(time, event) ~ x, ~ z, d
  )
  expect_equal(out$stan_data$N, 3L)
  expect_equal(out$stan_data$M, 3L)
  expect_equal(out$stan_data$K, 2L)  # intercept + x
  expect_equal(out$stan_data$Q, 2L)  # intercept + z
  expect_equal(out$stan_data$t0, rep(0, 3))
  expect_equal(sum(out$stan_data$spell_event), 2L)
})

test_that("group assignment creates correct J", {
  skip_if_not_installed("survival")
  d <- data.frame(time = 1:6, event = c(1,0,1,0,1,0),
                  x = rnorm(6), grp = rep(1:3, each = 2))
  out <- stanSPsurv:::prepare_spsurv_data(
    survival::Surv(time, event) ~ x, ~ 1, d, group = d$grp
  )
  expect_equal(out$stan_data$J, 3L)
  expect_equal(length(out$group_levels), 3L)
})

test_that("adjacency edges are extracted correctly", {
  A <- matrix(c(0,1,0, 1,0,1, 0,1,0), 3, 3)
  edges <- stanSPsurv:::adjacency_edges(A)
  expect_equal(edges$node1, c(1L, 2L))
  expect_equal(edges$node2, c(2L, 3L))
})

test_that("NULL adjacency returns empty edges", {
  edges <- stanSPsurv:::adjacency_edges(NULL)
  expect_equal(length(edges$node1), 0L)
})

test_that("invalid time intervals are rejected", {
  skip_if_not_installed("survival")
  d <- data.frame(time = c(-1, 2), event = c(0, 1), x = c(0, 1))
  expect_error(stanSPsurv:::prepare_spsurv_data(
    survival::Surv(time, event) ~ x, ~ 1, d
  ))
})

test_that("rank-deficient design matrix is rejected", {
  skip_if_not_installed("survival")
  d <- data.frame(time = c(1, 2, 3), event = c(0, 1, 0),
                  x1 = c(1, 2, 3), x2 = c(2, 4, 6))  # x2 = 2*x1
  expect_error(stanSPsurv:::prepare_spsurv_data(
    survival::Surv(time, event) ~ x1 + x2, ~ 1, d
  ), "full rank")
})
