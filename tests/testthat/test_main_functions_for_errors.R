context("Test that the main functions run well")


test_that("MixNRMI1 runs ok", {
  expect_output(str(MixNRMI1(rnorm(20), Nit = 50, extras = TRUE)), "List")
})
test_that("MixNRMI2 runs ok", {
  expect_output(str(MixNRMI2(rnorm(20), Nit = 50, extras = TRUE)), "List")
})
test_that("MixNRMI1cens runs ok", {
  expect_output(str(MixNRMI1cens(salinity$left, salinity$left, Nit = 50, extras = TRUE)), "List")
})
test_that("MixNRMI2cens runs ok", {
  expect_output(str(MixNRMI2cens(salinity$left, salinity$left, Nit = 50, extras = TRUE)), "List")
})