test_that("mixture density functions work", {
  expect_equal(
    pmix_vec_loop(
      qs = c(0.1, 0.2),
      locations_list = list(c(0, 0)),
      scales_list = list(c(1, 1)),
      weights_list = list(c(0.5, 0.5)),
      distr.k = 1
    ),
    pnorm(
      q = c(0.1, 0.2),
      mean = 0,
      sd = 1
    )
  )
  expect_equal(
    dmix_vec_loop(
      xs = c(0.1, 0.2),
      locations_list = list(c(0, 0)),
      scales_list = list(c(1, 1)),
      weights_list = list(c(0.5, 0.5)),
      distr.k = 1
    ),
    dnorm(
      x = c(0.1, 0.2),
      mean = 0,
      sd = 1
    )
  )
  expect_equal(
    qmix_vec_loop(
      ps = c(0.1, 0.2),
      locations_list = list(c(0, 0)),
      scales_list = list(c(1, 1)),
      weights_list = list(c(0.5, 0.5)),
      distr.k = 1
    ),
    qnorm(
      p = c(0.1, 0.2),
      mean = 0,
      sd = 1
    ),
    tolerance = .Machine$double.eps^0.25
  )
})
test_that("mixture censored density functions work", {
  expect_equal(
    dmixcens(
      xlefts = c(0.1, 0.2),
      xrights = c(0.1, 0.2),
      c_code_filters = list("1" = c(1, 2)),
      locations = c(0, 0),
      scales = c(1, 1),
      weights = c(0.5, 0.5),
      distr.k = 1
    ),
    dnorm(
      x = c(0.1, 0.2),
      mean = 0,
      sd = 1
    )
  )
})
test_that("numerical inversion quantile function do not produce errors", {
  expect_output(str(qmix(
    ps = c(0.1, 0.2),
    locations = c(0, 0),
    scales = c(1, 1),
    weights = c(0.5, 0.5),
    distr.k = 1
  )), "num [1:2]", fixed = TRUE)

  expect_output(str(qmix(
    ps = c(0.1, 0.2),
    locations = c(0.5, 0.5),
    scales = c(1, 1),
    weights = c(0.5, 0.5),
    distr.k = 2
  )), "num [1:2]", fixed = TRUE)

  expect_output(str(qmix(
    ps = c(0.1, 0.2),
    locations = c(0.5, 0.5),
    scales = c(0.4, 0.4),
    weights = c(0.5, 0.5),
    distr.k = 3
  )), "num [1:2]", fixed = TRUE)

  expect_output(str(qmix(
    ps = c(0.1, 0.2),
    locations = c(0, 0),
    scales = c(1, 1),
    weights = c(0.5, 0.5),
    distr.k = 4
  )), "num [1:2]", fixed = TRUE)

  expect_output(str(qmix(
    ps = c(0.1, 0.2),
    locations = c(0, 0),
    scales = c(1, 1),
    weights = c(0.5, 0.5),
    distr.k = 5
  )), "num [1:2]", fixed = TRUE)

  expect_output(str(qmix(
    ps = c(0.1, 0.9),
    locations = c(4.664140, 4.740357, 5.911687, 4.661989, 2.783218, 6.985979, 3.907708, 4.569024, 5.232156, 6.999026, 4.266442, 6.683336, 3.254375, 4.259333, 6.537262, 5.853631, 5.085109, 6.096547, 6.822146, 5.813071, 4.990496, 4.528902, 4.497958, 4.931694, 6.441307, 6.219373),
    scales = rep(0.2832576, 26),
    weights = c(
      0.0134999915, 0.0062952205, 0.0042048405, 0.0022923494, 0.0018878532, 0.0018635573, 0.0014329597, 0.0012816347, 0.0007893726, 0.0004881692,
      0.0004370247, 0.0002897057, 0.0004642882, 0.5012475440, 0.0376617422, 0.0148668278, 0.0749476782, 0.1437713219, 0.1342706424, 0.0185026008,
      0.0241549353, 0.0075760738, 0.0033351309, 0.0040153086, 0.0001638434, 0.0002593836
    ),
    distr.k = 1
  )), "num [1:2]", fixed = TRUE)
})

test_that("The vectorised mixture pdf calculation coincides with pmixnorm", {
  pmixnorm_vec_loop <- function(xs, means_list, sigmas_list, weights_list) {
    res <- 0.0 * xs
    nit <- length(means_list)
    for (it in 1:nit) {
      for (cmp in seq_along(means_list[[it]])) {
        res <- res + weights_list[[it]][cmp] * pnorm(q = xs, mean = means_list[[it]][cmp], sd = sigmas_list[[it]][cmp])
      }
    }
    return(res / nit)
  }
  data(acidity)
  xs <- seq(-5, 5, length.out = 100)
  outttest <- MixNRMI1(acidity, Nit = 50, extras = TRUE)
  outttest$sigmas_filled <- fill_sigmas(outttest)
  ref <- pmixnorm_vec_loop(xs = xs, means_list = outttest$means, sigmas_list = outttest$sigmas_filled, weights_list = outttest$weights)
  res <- pmix_vec_loop(qs = xs, locations_list = outttest$means, scales_list = outttest$sigmas_filled, weights_list = outttest$weights, distr.k = 1)
  expect_equal(res, ref)
})
