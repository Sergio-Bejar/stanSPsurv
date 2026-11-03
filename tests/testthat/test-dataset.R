test_that("conflict_spells dataset is available and well-formed", {
  data(conflict_spells, package = "stanSPsurv")
  expect_s3_class(conflict_spells, "data.frame")
  expect_equal(nrow(conflict_spells), 300L)
  expect_true(all(c("time", "event", "victory", "gdp", "peacekeeping",
                     "country", "region", "country_id") %in% names(conflict_spells)))
  expect_true(all(conflict_spells$time > 0))
  expect_true(all(conflict_spells$event %in% c(0L, 1L)))
  expect_equal(length(unique(conflict_spells$country)), 20L)
  expect_equal(length(unique(conflict_spells$region)), 4L)
})
