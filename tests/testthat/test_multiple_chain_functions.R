test_that("multMixNRMI1 runs ok", {
  expect_output(str(multMixNRMI1(rnorm(20), Nit = 50, extras = TRUE)), "List")
})
test_that("multMixNRMI2 runs ok", {
  expect_output(str(multMixNRMI2(rnorm(20), Nit = 50, extras = TRUE)), "List")
})
test_that("multMixNRMI1cens runs ok", {
  data(salinity)
  expect_output(str(multMixNRMI1cens(salinity$left, salinity$right, Nit = 50, extras = TRUE)), "List")
})
test_that("multMixNRMI2cens runs ok", {
  data(salinity)
  expect_output(str(multMixNRMI2cens(salinity$left, salinity$right, Nit = 50, extras = TRUE)), "List")
})