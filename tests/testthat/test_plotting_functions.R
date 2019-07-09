test_that("Plotting function for non censored data do not produce errors", {
  pdf(file = tempfile(pattern = "file", tmpdir = tempdir(), fileext = ""))
  data(acidity)
  outttest <- MixNRMI1(acidity, Nit = 5, extras = TRUE)
  p = plotCDF_noncensored(outttest)
  expect_output(str(p), "List of 9")
  p = plotPDF_noncensored(outttest)
  expect_output(str(p), "List of 9")
  p = qq_plot_noncensored(outttest)
  expect_output(str(p), "List of 9")
  p = plotGOF(outttest)
  expect_output(str(p), "gtable")

  outttest2 <- MixNRMI2(acidity, Nit = 5, extras = TRUE)
  p = plotCDF_noncensored(outttest2)
  expect_output(str(p), "List of 9")
  p = plotPDF_noncensored(outttest2)
  expect_output(str(p), "List of 9")
  p = qq_plot_noncensored(outttest2)
  expect_output(str(p), "List of 9")
  p = plotGOF(outttest2)
  expect_output(str(p), "gtable")

  data(enzyme)
  outttest3 <- MixNRMI2(enzyme, Alpha = 1, Kappa = 0.007, Gama = 0.5,
                                                  distr.k = 2, distr.py0 = 2,
                                                  distr.pz0 = 2, mu.pz0 = 1, sigma.pz0 = 1, Meps=0.005,
                                                  Nit = 50, Pbi = 0.2, extras = TRUE)
  p = plotCDF_noncensored(outttest3)
  expect_output(str(p), "List of 9")
  p = plotPDF_noncensored(outttest3)
  expect_output(str(p), "List of 9")
  p = plotGOF(outttest3)
  expect_output(str(p), "gtable")
  dev.off()
})

test_that("Plotting function for censored data do not produce errors", {
  pdf(file = tempfile(pattern = "file", tmpdir = tempdir(), fileext = ""))
  data(salinity)
  outttest <- MixNRMI1cens(salinity$left, salinity$right, Nit = 5, extras = TRUE)
  p = plotCDF_censored(outttest)
  expect_output(str(p), "List of 9")
  p = plotPDF_censored(outttest)
  expect_output(str(p), "List of 9")
  p = qq_plot_censored(outttest)
  expect_output(str(p), "List of 9")
  p = plotGOF(outttest)
  expect_output(str(p), "gtable")
  dev.off()
})

test_that("The vectorised mixture pdf calculation coincides with pmixnorm", {
  pmixnorm_vec_loop = function(xs, means_list, sigmas_list, weights_list){
    res = 0.0*xs
    nit = length(means_list)
    for (it in 1:nit){
      for (cmp in seq_along(means_list[[it]]))
        res = res + weights_list[[it]][cmp] * pnorm(q = xs, mean = means_list[[it]][cmp], sd = sigmas_list[[it]][cmp])
    }
    return(res/nit)
  }
  data(acidity)
  xs = seq(-5,5,length.out = 100)
  outttest <- MixNRMI1(acidity, Nit = 50, extras = TRUE)
  outttest$sigmas_filled = fill_sigmas(outttest)
  ref = pmixnorm_vec_loop(xs = xs, means_list = outttest$means, sigmas_list = outttest$sigmas_filled, weights_list = outttest$weights)
  res = pmix_vec_loop(qs = xs, locations_list = outttest$means, scales_list = outttest$sigmas_filled, weights_list = outttest$weights, distr.k = 1)
  expect_equal(res, ref)
})
