test_that("The kernel name function works", {
  expect_equal(as.character(BNPdensity:::give_kernel_name(4)), "Double Exponential")
})
