test_that("mixture density functions work", {
  expect_equal(
    pmix_vec_loop(
      xs = c(0.1, 0.2),
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
