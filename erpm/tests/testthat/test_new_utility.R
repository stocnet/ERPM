### Tests Utility functions

# Size statistics for the 3
test_that("Statistics size", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  expect_equal(stat_size(p1,'avg'), 1)
  expect_equal(stat_size(p2,'avg'), 6)
  expect_equal(stat_size(p3,'avg'), 1.5)

  expect_equal(stat_size(p1,'sd'), 0)
  expect_true(is.na(stat_size(p2,'sd')))
  expect_equal(round(stat_size(p3,'sd'),3), 0.577)
})


## Proportion of isolates ---------------------------------

test_that("Proportion of isolates", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  expect_equal(prop_isolate(p1), 1)
  expect_equal(prop_isolate(p2), 0)
  expect_equal(prop_isolate(p3), 1/3)

})


## Range ---------------------------------

# All the same

test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  at <- rep(1,6)


  expect_equal(range(p1,at,'sum'), 0)
  expect_equal(range(p2,at,'sum'), 0)
  expect_equal(range(p3,at,'sum'),0)

})

# Mixed attribute

test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  at <- c(2,3,2,1,3,4)


  expect_equal(range(p1,at,'sum'), 0)
  expect_equal(range(p2,at,'sum'), 3)
  expect_equal(range(p3,at,'sum'),3)

})


# Averager per group

# All the same
test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  at <- rep(1,6)


  expect_equal(range(p1,at,'avg'), 0)
  expect_equal(range(p2,at,'avg'), 0)
  expect_equal(range(p3,at,'avg'),0)

})

# Mixed attribute

test_that("Range", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  at <- c(2,3,2,1,3,4)


  expect_equal(range(p1,at,'avg'), 0)
  expect_equal(range(p2,at,'avg'), 3)
  expect_equal(range(p3,at,'avg'),3/4)

})


## ICC ---------------------------------

# All the same

test_that("ICC", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  at <- rep(1,6)


  expect_true(is.na(icc(p1,at))) # NAN problem
  expect_equal(icc(p2,at), 1)
  expect_equal(icc(p3,at),1)

})


# Mixed attribute

test_that("icc", {
  p1 <- 1:6
  p2 <- rep(1,6)
  p3 <- c(1,2,2,3,3,4)

  at <- c(2,3,2,1,3,4)


  expect_true(is.na(icc(p1,at)))
  expect_true(is.na(icc(p2,at)))
  expect_equal(round(icc(p3,at),2),0.42)

})

