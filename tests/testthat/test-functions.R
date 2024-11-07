test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("sample_proportion_prob works correctly", {
  result <- sample_proportion_prob(0.5, 100, 0.55)
  expect_type(result, "double")
})
