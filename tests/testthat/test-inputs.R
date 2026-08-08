test_that("adjacency validation rejects malformed matrices", {
  expect_error(stanSPsurv:::validate_adjacency(diag(2), 2), "zero diagonal")
  expect_silent(stanSPsurv:::validate_adjacency(matrix(c(0,1,1,0), 2), 2))
})

