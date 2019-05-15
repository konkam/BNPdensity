context("Test the plotting functions")

#' data(acidity)
#' out <- MixNRMI1(acidity, extras = TRUE)
#' plotCDF_noncensored(out)
#'

test_that("The CDF plot does not produce ane error", {
  data(acidity)
  out <- MixNRMI1(acidity, Nit = 50, extras = TRUE)
  expect_output(str(plotCDF_noncensored(out)), "List of 9")
})