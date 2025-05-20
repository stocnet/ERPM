### Tests Utility functions

# Size statistics for the 3
test_that("Statistics size", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  expect_equal(group_size(p1, "avg"), 1)
  expect_equal(group_size(p2, "avg"), 6)
  expect_equal(group_size(p3, "avg"), 1.5)

  expect_equal(group_size(p1, "sd"), 0)
  expect_true(is.na(group_size(p2, "sd")))
  expect_equal(round(group_size(p3, "sd"), 3), 0.577)
})


## Proportion of isolates ---------------------------------

test_that("Proportion of isolates", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  expect_equal(proportion_isolate(p1), 1)
  expect_equal(proportion_isolate(p2), 0)
  expect_equal(proportion_isolate(p3), 1 / 3)

})


## Range ---------------------------------

# All the same

test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  at <- rep(1, 6)


  expect_equal(range_attribute(p1, at, "sum_pergroup"), 0)
  expect_equal(range_attribute(p2, at, "sum_pergroup"), 0)
  expect_equal(range_attribute(p3, at, "sum_pergroup"), 0)

})

# Mixed attribute

test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  at <- c(2, 3, 2, 1, 3, 4)


  expect_equal(range_attribute(p1, at, "sum_pergroup"), 0)
  expect_equal(range_attribute(p2, at, "sum_pergroup"), 3)
  expect_equal(range_attribute(p3, at, "sum_pergroup"), 3)

})


# Average per group

# All the same
test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  at <- rep(1, 6)


  expect_equal(range_attribute(p1, at, "avg_pergroup"), 0)
  expect_equal(range_attribute(p2, at, "avg_pergroup"), 0)
  expect_equal(range_attribute(p3, at, "avg_pergroup"), 0)

})

# Mixed attribute

test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  at <- c(2, 3, 2, 1, 3, 4)


  expect_equal(range_attribute(p1, at, "avg_pergroup"), 0)
  expect_equal(range_attribute(p2, at, "avg_pergroup"), 3)
  expect_equal(range_attribute(p3, at, "avg_pergroup"), 3 / 4)

})


## ICC ---------------------------------

# All the same attribute

test_that("ICC", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  at <- rep(1, 6)


  expect_true(is.na(icc(p1, at))) # NAN problem
  expect_true(is.na(icc(p2, at))) # NAN problem
  expect_true(is.na(icc(p3, at))) # NAN problem

})


# Mixed attribute

test_that("icc", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  at <- c(2, 3, 2, 1, 3, 4)

  expect_true(is.na(icc(p1, at)))
  expect_true(is.na(icc(p2, at)))
  expect_equal(icc(p3, at), 11 / 26)

})


# Different attribute

test_that("icc", {
  p1 <- 1:6
  p2 <- rep(1, 6)
  p3 <- c(1, 2, 2, 3, 3, 4)

  at <- 1:6

  expect_true(is.na(icc(p1, at)))
  expect_true(is.na(icc(p2, at)))
  expect_equal(icc(p3, at), 29 / 32)

})
