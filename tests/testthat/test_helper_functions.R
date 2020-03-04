test_that("Correct testing for censored data", {
  data(salinity)
  expect_true(is_censored(salinity))
  data(acidity)
  expect_false(is_censored(acidity))
})

test_that("Correct extraction of CPOs", {
  data(acidity)
  out <- MixNRMI1(acidity, Nit = 50)
  cpos = cpo(out)
  expect_equal(cpos, out$cpo)
  out <- MixNRMI2(acidity, Nit = 50)
  cpos = cpo(out)
  expect_equal(cpos, out$cpo)
  out <- multMixNRMI1(acidity, parallel = TRUE, Nit = 10, ncores = 2)
  cpos = cpo(out)
  expect_equal(length(cpos), length(acidity))
  })
