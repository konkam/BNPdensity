testthat::test_that("Conversion to mcmc object works", {
  out <- BNPdensity::convert_to_mcmc(BNPdensity::multMixNRMI1(rnorm(20), Nit = 50, extras = TRUE))
  testthat::expect_is(out, "mcmc")
  BNPdensity::data("salinity")
  out <- BNPdensity::convert_to_mcmc(BNPdensity::multMixNRMI1cens(salinity$left, salinity$right, Nit = 50, extras = TRUE))
  testthat::expect_is(out, "mcmc")
})
