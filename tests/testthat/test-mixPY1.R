test_that("MixPY1 runs ok", {
  expect_output(str(MixPY1(x = rnorm(20))), "List")
})
