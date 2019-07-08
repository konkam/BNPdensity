test_that("Correct testing for censored data", {
  data(salinity)
  expect_true(is_censored(salinity))
  data(acidity)
  expect_false(is_censored(acidity))
})
