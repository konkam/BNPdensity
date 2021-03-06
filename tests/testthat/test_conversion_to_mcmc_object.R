test_that("Conversion to mcmc object works", {
  out <- convert_to_mcmc(multMixNRMI1(rnorm(20), Nit = 50, extras = TRUE))
  expect_is(out, "mcmc.list")
  data("salinity")
  out <- convert_to_mcmc(multMixNRMI1cens(salinity$left, salinity$right, Nit = 50, extras = TRUE))
  expect_is(out, "mcmc.list")
  library(coda)
  out <- as.mcmc.list(as.mcmc(multMixNRMI1cens(salinity$left, salinity$right, Nit = 50, extras = TRUE)))
  expect_is(out, "mcmc.list")
})

test_that("The coda methods work", {
  data(acidity)
  out <- convert_to_mcmc(multMixNRMI1(acidity, Nit = 50, extras = TRUE))
  expect_silent(plot(out))
  expect_silent(acfplot(out))
  expect_silent(gelman.diag(out))
  out <- as.mcmc.list(as.mcmc(multMixNRMI1cens(salinity$left, salinity$right, Nit = 50, extras = TRUE)))
  expect_silent(plot(out))
  expect_silent(acfplot(out))
  expect_silent(gelman.diag(out))
})
