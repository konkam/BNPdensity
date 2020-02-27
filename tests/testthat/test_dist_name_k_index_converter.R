test_that("Distribution name conversion function works", {
  expect_equal(as.numeric(dist_name_k_index_converter('norm')), 1)
  expect_error(as.numeric(dist_name_k_index_converter(6)))
  expect_equal(as.numeric(process_dist_name('norm')), 1)
  expect_equal(as.numeric(process_dist_name(1)), 1)
})
