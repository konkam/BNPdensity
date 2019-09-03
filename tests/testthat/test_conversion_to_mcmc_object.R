test_that("Conversion to mcmc object works", {
  out <- convert_to_mcmc(multMixNRMI1(rnorm(20), Nit = 50, extras = TRUE))
  expect_is(out, "mcmc")
  data("salinity")
  out <- convert_to_mcmc(multMixNRMI1cens(salinity$left, salinity$right, Nit = 50, extras = TRUE))
  expect_is(out, "mcmc")
  library(coda)
  out <- as.mcmc(multMixNRMI1cens(salinity$left, salinity$right, Nit = 50, extras = TRUE))
  expect_is(out, "mcmc.list")
})

